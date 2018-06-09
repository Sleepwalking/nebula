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

function y = preprocess(x, dither_level = 0.05, dc_cutoff = 50 / 4000)
  xsqr_intg = cumsum(x .^ 2);
  xrms = sqrt((xsqr_intg(257:end) - xsqr_intg(1:end - 256)) / 256);
  xrms = [ones(128, 1) * xrms(1); xrms; ones(128, 1) * xrms(1)];
  peak = max(xrms);
  thrd = peak * dither_level;
  x = dcnotch(x, dc_cutoff);
  y = x + (xrms < thrd) .* randn(size(x)) .* (thrd - xrms);
end

