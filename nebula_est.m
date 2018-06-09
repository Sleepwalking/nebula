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

function [f0, v, pv, lmap] = nebula_est(model, x, fs, thop = 0.005)
  if fs ~= 8000
    x = downsample(x(:, 1), fs, 8000);
  end
  x = preprocess(x);
  fs = 8000;
  
  nch = length(model.G);
  taxis = (0:thop:(length(x) - 1) / fs)';
  naxis = ceil(taxis * fs + eps);
  fc = logspace(log10(40), log10(1000), nch)' / fs;
  fres = fc;
  
  [SNR1 IF1 SNR2 IF2 SNR0] = estimate_ifsnr( ...
    x, naxis, fc, fres, 'nuttall98', 'hanning');

  eps = 1e-10;
  SNR0 = mag2db(max(SNR0, eps));
  SNR1 = mag2db(max(SNR1, eps));
  SNR2 = mag2db(max(SNR2, eps));
  IF1 = IF1 * fs;
  IF2 = IF2 * fs;

  [Lmap f] = make_likelihood_map(model, SNR1, IF1, SNR2, IF2, SNR0);
  lmap = postprocess(exp(Lmap), 1.5 ./ f / fs / thop);

  LF0 = repmat(1:columns(Lmap), rows(Lmap), 1);
  ptrans = normpdf((0:columns(Lmap-1)), 0, 2);
  [s Ltotal] = vitfilt1(LF0, lmap, ptrans / sum(ptrans));

  [q, pv] = detect_voicing_status(log(lmap), s);
  v = 2 - q;
  f0 = refine_f0(log(lmap), s, f);
  if nargout() == 1
    f0 = f0 .* v;
  end
end

function [Lmap, f] = make_likelihood_map(model, SNR1, IF1, SNR2, IF2, SNR0)
  fs = 8000;
  nch = length(model.G);
  nf = length(model.Lcal);
  f = logspace(log10(40), log10(1000), nf)' / fs;
  Lmap = zeros(rows(SNR1), nf) - 100;
  Lsub = {};
  for i = 1:nch
    Lsub{i} = zeros(rows(SNR1), nf);
  end
  for i = 1:rows(SNR1)
    Li = zeros(nch, columns(Lmap));
    for j = 1:nch
      L = condgmm(model.G{j}, f * fs, ...
        [SNR1(i, j), SNR0(i, j), SNR2(i, j), IF1(i, j), IF2(i, j)]', 2);
      L -= mean(L);
      Lsub{j}(i, :) = L;
      Li(j, :) = L;
    end
    Lmap(i, :) = (mean(Li) - mean(model.Lcal));
    Lmap(i, :) -= log(sum(exp(Lmap(i, :))));
  end
end

function [q pvoiced] = detect_voicing_status(Lmap, s)
  f0p = zeros(rows(Lmap), 1);
  for i = 1:rows(Lmap)
    f0p(i) = Lmap(i, s(i));
  end
  umean = -4.78;
  ustd = 0.12;

  Lfunc = @(x) -mean( ...
    log(normpdf(f0p, x(1), x(2)) + normpdf(f0p, umean, ustd)));
  param = fminunc(Lfunc, [umean(1) + 4, 1]);
  pvoiced = normpdf(f0p, param(1), param(2)) ./ ...
    (normpdf(f0p, param(1), param(2)) + normpdf(f0p, umean, ustd));
  [q, ptotal] = binvitsearch(pvoiced, 0.01);
end

function f0 = refine_f0(Lmap, s, f)
  nf = length(f);
  Lidx = 1:0.05:nf;
  Linterp = interp1(1:nf, Lmap', Lidx, 'spline')';
  sinterp = interp1(Lidx, 1:length(Lidx), s);

  f0 = zeros(rows(s), 1);
  for i = 1:rows(s)
    iidx = max(1, sinterp(i)-10:sinterp(i)+10);
    [~, is] = max(Linterp(i, iidx));
    f0(i) = iidx(is);
  end
  f0 = interp1(1:nf, f, f0 / 20) * 8000;
end

