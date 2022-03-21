function [Tf,Tfm] = flipheatfix(t,h0,h1,Nm,Nz)
%FLIPHEATFIX   The fixed point of the flip-heat operation for fixed time.
%   TF = FLIPHEATFIX(T,H0,H1,NM,NZ) returns the the temperature profile TF
%   that is the fixed point of the flip-heat operator with flip duration T.
%   See <strong>heateigfun</strong> for a description of H0,H1,NM,NZ.
%
%   TF = FLIPHEATFIX(T,H0,H1,IFT,MU) uses the matrix IFT and vector MU
%   returned by <strong>heateigfun</strong> rather than computing them from
%   scratch.
%
%   [TF,TFM] = FLIPHEATFIX(...) also returns the Fourier components of TF.
%
%   See also FLIPOP, HEATEIGFUN.

%
% This file is part of cookflip_code
%
% Copyright (c) 2022 Jean-Luc Thiffeault <jeanluc@math.wisc.edu>
%
% See the file LICENSE for copying permission.
%

if nargin < 1 || isempty(t), t = .02; end
% Use "cooking" values if h0, h1 not given.
if nargin < 2 || isempty(h0), h0 = 21.6; end
if nargin < 3 || isempty(h1), h1 = 1.44; end
if nargin < 4 || isempty(Nm), Nm = 31; end
if nargin < 5 || isempty(Nz), Nz = 1001; end

%% Lots of duplicate code from cooktime.

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

% Compute by integrating over several periods.  More accurate?

% Equilibrium profile and its Fourier coefficients.
Teq = heatsteady(h0,h1,Nz); Teqm = heatsteady(h0,h1,mu.');

% The "flipping" matrix: project flipped modes onto unflipped modes.
f = flipop(IFT);

D = diag(exp(-mu.^2*t));
F = eye(Nm) - D*f;

if false
  % This is a direct but bad way to compute the fixed profile.  Lots of
  % aliasing.
  Tf = IFT.' * (F \ (eye(Nm) - D)*Teqm);
  return
end

% The flipped equilibrium profile.
B = Teqm - f*Teqm;

% The symmetrized decay matrix.
%scriptF = f'*D.^2*f
%-1/log(norm(scriptF))  % how many intervals to decrease by 1/e.

F = eye(Nm) - D*f;

thetafpm = -F\D*B;
thetafp = IFT'*thetafpm;

Tf = Teq + thetafp;
Tfm = Teqm + thetafpm;
