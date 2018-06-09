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

function y = dcnotch(x, cutoff)
  a1 = - 2.0 * cos(pi * cutoff);
  a0 = 8.0 * cos(pi * cutoff) - 7.0;
  r = (-a1 - sqrt(a1 ^ 2 - 4.0 * a0)) / 2.0;
  a = [1.0, -r];
  b = [1.0, -1.0];
  y = filtfilt(b, a, x);
end

