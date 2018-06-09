% Nebula
% ===
% Copyright (c) 2018 Kanru Hua. All rights reserved.

% This file is part of Nebula.

% Nebula is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Nebula is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with Nebula. If not, see <http://www.gnu.org/licenses/>.

wtype = 'nuttall98';
stype = 'hanning';

nch = 36;
data_dir = 'data';
model_dir = 'model';

fs = 8000;
nch = 36;
nf = 128;
fc = logspace(log10(40), log10(1000), nch)' / fs;
f = logspace(log10(40), log10(1000), nf)' / fs;
fres = fc;

[x naxis] = make_random_signal(f0, 100000, -100);
[SNR1 IF1 SNR2 IF2 SNR0] = estimate_ifsnr(x, naxis, fc, fres, wtype, stype);

eps = 1e-10;
SNR0 = mag2db(max(SNR0, eps));
SNR1 = mag2db(max(SNR1, eps));
SNR2 = mag2db(max(SNR2, eps));
IF1 = IF1 * fs;
IF2 = IF2 * fs;

Lcal = zeros(nch, nf);
for j = 1:nch
  system(sprintf( ...
    'gmm -m 16 -l 6 -v 0.001 -w 0.001 -f %s/samples-%d.f > %s/gmm16-%d.f', ...
    data_dir, j, model_dir, j));
  G = load_gmm(sprintf('%s/gmm16-%d.f', model_dir, j), 16, 6);
  Lj = zeros(rows(SNR0), nf);
  for i = 1:rows(SNR0)
    L = condgmm(G, f * fs, ...
      [SNR1(i, j), SNR0(i, j), SNR2(i, j), IF1(i, j), IF2(i, j)]', 2);
    L -= mean(L);
    Lj(i, :) = L;
  end
  Lcal(j, :) = mean(Lj);
end

save(sprintf('%s/Lcal', model_dir), 'Lcal');

