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

function [x naxis] = make_random_signal(f0, n, at)
  xn = (1:n)';
  x = zeros(size(xn));
  nhar = floor(0.5 / f0);
  for i = 1:nhar
    p = rand * 2 * pi;
    h = 10 ^ (rand * 1.0 - 0.5);
    if i == 1
      h = 1;
    end
    x += h * sin(2 * pi * f0 * i * xn + p);
  end
  x += randn(size(x)) * at;
  naxis = (1:512:length(x))';
end

