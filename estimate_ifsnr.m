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

% estimate_ifsnr.m
% Instantaneous Frequency (IF) and Signal-to-Noise Ratio (SNR) feature extractor.
% Inputs
%   x       column vector       input signal
%   naxis   column vector       analysis time indices
%   fc      column vector       central frequencies (as ratio of sampling rate)
%   fres    column vector       band width (as ratio of sampling rate)
%   wtype   string              type of analysis window (see fucntion getwcoef)
%   stype   string              type of smoothing window (see fucntion getwcoef)
% Outputs
%   SNR1    matrix              SNR feature
%   SNR2    matrix              SNR feature at double central frequency
%   SNR0    matrix              SNR feature at half central frequency
%   IF1     matrix              IF feature
%   IF2     matrix              IF feature at double central frequency
function [SNR1 IF1 SNR2 IF2 SNR0] = estimate_ifsnr( ...
    x, naxis, fc, fres, wtype = 'nuttall83', stype = 'nuttall83')
  nout = nargout();
  nch  = length(fc);
  SNR1 = zeros(length(naxis), nch);
  SNR2 = zeros(length(naxis), nch);
  IF1  = zeros(length(naxis), nch);
  IF2  = zeros(length(naxis), nch);
  
  aw   = getwcoef(wtype);
  aws  = getwcoef(stype);
  Y    = zeros(length(x), nch);

  for i = 1:nch
    nw   = ceil(4 / fres(i));
    nw_2 = floor(nw / 2);
    omegaw = 2 * pi * fres(i) / 4;
    w = zeros(nw, 1);
    ws = zeros(nw, 1);
    for j = 1:4
      jw  = cos((j - 1) * omegaw * ((0:nw-1)' - nw_2));
      w  += aw (j) * jw;
      ws += aws(j) * jw;
    end
    w  /= sum(w);
    ws /= sum(ws);

    omegah = 2 * pi * fc(i);
    h  = w .* exp(J * omegah * ((0:nw-1)' - nw_2));
    h2 = w .* exp(J * omegah * ((0:nw-1)' - nw_2) * 2);

    SNR1(:, i)   = x2as(x, h , ws)(naxis);
    if nout > 2
      SNR2(:, i) = x2as(x, h2, ws)(naxis);
    end

    fi = ifdetect(x, naxis, fc(i), fres(i));
    fi = min(max(fi, fc(i) * 0.5), fc(i) * 1.5);
    IF1(:, i) = fi;

    if nout > 2
      fi = ifdetect(x, naxis, fc(i) * 2, fres(i));
      fi = min(max(fi, fc(i) * 1.5), fc(i) * 2.5);
      IF2(:, i) = fi;
    end
  end
  
  SNR0 = [SNR1(:, 8:14) SNR1(:, 1:end-7)];
end

function as = x2as(x, h, ws)
  nw   = length(h);
  nw_2 = floor(nw / 2);
  y1   = fftconv(x, h)(nw_2+1:end-nw+nw_2+1);
  y1 ./= abs(y1) + eps;
  y2   = fftconv(y1, h)(nw_2+1:end-nw+nw_2+1);
  r    = y1 - y2;
  a    = abs(r) .^ 2;
  as   = fftconv(a, ws)(nw_2+1:end-nw+nw_2+1);
  as(1:nw)       = as(nw);
  as(end-nw:end) = as(end-nw);
end

function a = getwcoef(wtype)
  if strcmp(wtype, 'nuttall83')
    a = [0.338946, 0.481973, 0.161054, 0.018027];
  elseif strcmp(wtype, 'nuttall98')
    a = [0.3635819, 0.4891775, 0.1365995, 0.0106411];
  elseif strcmp(wtype, 'nuttall93')
    a = [0.355768, 0.487396, 0.144232, 0.012604];
  elseif strcmp(wtype, 'nuttall64')
    a = [0.40897, 0.5, 0.09103, 0];
  elseif strcmp(wtype, 'hanning')
    a = [0.5, 0.5, 0, 0];
  elseif strcmp(wtype, 'hamming')
    a = [0.54, 0.46, 0, 0];
  elseif strcmp(wtype, 'blackman')
    a = [0.42, 0.5, 0.08, 0];
  elseif strcmp(wtype, 'blackman-harris')
    a = [0.4243801, 0.4973406, 0.0782793, 0];
  elseif strcmp(wtype, 'boxcar')
    a = [1, 0, 0, 0];
  end
end

