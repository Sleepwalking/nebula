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

pkg load control;
pkg load signal;

wtype = 'nuttall98';
stype = 'hanning';

fs = 8000;
nch = 36;
fc = logspace(log10(40), log10(1000), nch)' / fs;
fres = fc;

data = {};
num_samples = 0;

eps = 1e-10;
for i = 1:3600
  printf('%d/%d\n', i, 100); fflush(stdout);
  
  global_snr = db2mag(rand * 100 - 50);
  f0 = (40 + (1000 - 40) * rand) / fs;
  [x naxis] = make_random_signal(f0, 10000, global_snr);
  
  [SNR1 IF1 SNR2 IF2 SNR0] = estimate_ifsnr(x, naxis, fc, fres, wtype, stype);
  
  data{i}.SNR0 = mag2db(max(SNR0, eps));
  data{i}.SNR1 = mag2db(max(SNR1, eps));
  data{i}.SNR2 = mag2db(max(SNR2, eps));
  data{i}.IF1 = IF1 * fs;
  data{i}.IF2 = IF2 * fs;
  data{i}.f0 = f0 * fs;
  
  num_samples += rows(SNR0);
end

for i = 1:nch
  samples = zeros(num_samples, 6);
  idx = 1;
  for j = 1:length(data)
    num = rows(data{j}.SNR0);
    F0 = data{j}.f0 * ones(num, 1);
    samples(idx:idx + num - 1, :) = [ ...
      data{j}.SNR1(:, i) data{j}.SNR0(:, i) data{j}.SNR2(:, i) ...
      data{j}.IF1(:, i) data{j}.IF2(:, i) F0];
    idx += num;
  end
  fd = fopen(sprintf('data/samples-%d.f', i), 'wb');
  fwrite(fd, samples'(:), 'float');
  fclose(fd);
end

