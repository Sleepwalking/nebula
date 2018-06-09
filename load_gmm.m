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

% load_gmm
% Load model parameters from binary GMM file (in SPTK format)
% Inputs
%   filename    string          path to the GMM file
%   nmix        integer         number of mixture components
%   nord        integer         dimensionality
% Output
%   G           cell array      GMM model parameters
%   G{i}        structure       model parameters for the i-th mixture component
%   G{i}.mu     column vector   mean vector
%   G{i}.sigma  matrix          covariance matrix
%   G{i}.inv    matrix          inverse of the covariance matrix
%   G{i}.det    real            determinant
%   G{i}.inv11  matrix          inverse of the (1, 1) block cov. matrix
%   G{i}.det11  matrix          determinant of the (1, 1) block cov. matrix
%   G{i}.condsigma  real        variance of the conditional on the (2, 2) block
%   G{i}.s21s11i    row vector  product of the last row of the (2, 1) block cov.
%                               matrix with the inverse of the (1, 1) block cov.
%                               matrix
function G = load_gmm(filename, nmix, nord)
  G = {};

  fp = fopen(filename, 'rb');
  x = fread(fp, inf, 'float');
  fclose(fp);

  for i = 1:nmix
    G{i}.w = x(i);
  end

  x = x(nmix+1:end);
  for i = 1:nmix
    sigma = reshape(x(nord + 1:nord + nord * nord), nord, nord);
    G{i}.mu        = x(1:nord);
    G{i}.sigma     = sigma;
    G{i}.inv       = inv(sigma);
    G{i}.det       = det(sigma);
    G{i}.inv11     = inv(sigma(1:end - 1, 1:end - 1));
    G{i}.det11     = det(sigma(1:end - 1, 1:end - 1));
    G{i}.condsigma = sigma(end, end) - ...
      sigma(end, 1:end - 1) * G{i}.inv11 * sigma(1:end - 1, end);
    x = x(nord+nord*nord+1:end);
    G{i}.s21s11i   = sigma(end, 1:end - 1) * G{i}.inv11;
  end
end

