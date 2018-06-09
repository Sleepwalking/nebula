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

function [q, ptotal] = binvitsearch(pob, ptrans)
  a = zeros(length(pob), 2);
  bt = zeros(length(pob), 2);
  pstay = log(1 - ptrans);
  ptrans = log(ptrans);
  pob1 = log(pob);
  pob2 = log(1 - pob);
  a(1, 1) = pob1(1);
  a(1, 2) = pob2(1);
  for i = 2:length(pob)
    if a(i - 1, 2) + ptrans > a(i - 1, 1) + pstay
      a(i, 1) = a(i - 1, 2) + ptrans + pob1(i);
      bt(i, 1) = 2;
    else
      a(i, 1) = a(i - 1, 1) + pstay + pob1(i);
      bt(i, 1) = 1;
    end
    if a(i - 1, 1) + ptrans > a(i - 1, 2) + pstay
      a(i, 2) = a(i - 1, 1) + ptrans + pob2(i);
      bt(i, 2) = 1;
    else
      a(i, 2) = a(i - 1, 2) + pstay + pob2(i);
      bt(i, 2) = 2;
    end
  end
  [ptotal, last] = max(a(end, :));
  q = zeros(length(pob), 1);
  q(end) = last;
  for i = length(q)-1:-1:1
    q(i) = bt(i + 1, q(i + 1));
  end
end
