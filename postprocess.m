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

function Lsmooth = postprocess(L, navg)
  Lsmooth = L;
  for i = 1:columns(L)
    cs = cumsum(L(:, i)) / navg(i) / 2;
  
    idx = (0:rows(L)-1)';
    idx_h = idx + navg(i);
    idx_l = idx - navg(i);
    Lsmooth(:, i) = interp1(idx, cs, idx_h, 'extrap') - ...
                    interp1(idx, cs, idx_l, 'extrap');
  end
end

