function T = heat(t,T0,h0,h1,Nm,Nz)
%HEAT   Temperature profile for the heat equation.
%   T = HEAT(t,T0,H0,H1,NM,NZ) returns the the temperature profile T solving
%   the heat equation for a time t.  See <strong>heateigfun</strong> for a
%   description of H0,H1,NM,NZ.
%
%   TF = HEAT(t,T0,H0,H1,IFT,MU) uses the matrix IFT and vector MU returned
%   by <strong>heateigfun</strong> rather than computing them from scratch.
%
%   See also FLIPOP, HEATEIGFUN.

%
% This file is part of cookflip_code
%
% Copyright (c) 2022 Jean-Luc Thiffeault <jeanluc@math.wisc.edu>
%
% See the file LICENSE for copying permission.
%

if nargin < 1 || isempty(t), t = .1; end
if nargin < 2 || isempty(T0), T0 = 0; end
% Use "cooking" values if h0, h1 not given.
if nargin < 3 || isempty(h0), h0 = 21.6; end
if nargin < 4 || isempty(h1), h1 = 1.44; end
if nargin < 5 || isempty(Nm), Nm = 31; end
if nargin < 6 || isempty(Nz), Nz = 1001; end

if isscalar(Nm)
  % Compute IFT matrix.
  [IFT,mu] = heateigfun(h0,h1,Nm,Nz);
else
  % IFT and mu as arguments.
  IFT = Nm; mu = Nz;
  Nm = size(IFT,1); Nz = size(IFT,2);
end

m = 1:Nm;

if ~isinf(Nz)
  z = linspace(0,1,Nz);
  intcols = @(f) trapz(z,f,2);     % integral of columns
else
  intcols = @(f) integral(f);      % integral of columns
end

% Equilibrium profile and its Fourier coefficients.
Teq = heatsteady(h0,h1,Nz); Teqm = heatsteady(h0,h1,mu.');
if T0 ~= 0
  if ~isinf(Nz)
    T0m = zeros(Nm,1);
    for i = 1:Nm
      T0m(i) = intcols(IFT(i,:).*T0.');
    end
  else
    T0m = IFT*T0;
  end
else
  T0m = 0;
end

T = Teq + IFT.'*diag(exp(-mu.^2*t))*(T0m - Teqm);
