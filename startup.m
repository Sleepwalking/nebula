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

addpath('./mex');
addpath('./train');

if exist('condgmm') ~= 3
  system(sprintf('mkoctfile --mex mex/condgmm.c -o mex/condgmm.mex'));
end

if exist('ifdetect') ~= 3
  system(sprintf('mkoctfile --mex mex/ifdetect.c -o mex/ifdetect.mex'));
end

if exist('vitfilt1') ~= 3
  system(sprintf('mkoctfile --mex mex/vitfilt1.c external/libgvps/build/libgvps.a -o mex/vitfilt1.mex'));
end

