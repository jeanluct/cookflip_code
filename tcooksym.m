function [tc,tc1] = tcooksym(Tcook,h0,Nm,Nz)
%TCOOKSYM   Time to cook food for equal flips and symmetric boundary conditions.
%   TC = TCOOKSYM(COOKTEMP,H,NM,NZ) returns the time TC needed to cook a
%   slab of food of unit thickness, with equal-timed flips, when the boundary
%   conditions are symmeytric (H0=H1=H).  TC=NaN if this is not possible
%   (i.e., COOKTEMP >= .5, the midpoint temperature).  See <strong>heateigfun</strong> for a
%   description of H0,H1,NM,NZ.
%
%   TC = TCOOKSYM(COOKTEMP,H,IFT,MU) uses the matrix IFT and vector MU
%   returned by <strong>heateigfun</strong> rather than computing them from scratch.
%
%   [TC,TC1] = TCOOKSYM(...) also returns TC1, the 1-mode approximation
%   to the cooking time.
%
%   See also HEATEIGFUN, TEMP, TCOOKTHRU.

%
% This file is part of cookflip_code
%
% Copyright (c) 2022 Jean-Luc Thiffeault <jeanluc@math.wisc.edu>
%
% See the file LICENSE for copying permission.
%

%% This contains quite a bit of duplicate code with tcookthru.

if nargin < 1 || isempty(Tcook), Tcook = .257; end
if nargin < 2 || isempty(h0), h0 = 21.6; end
if nargin < 3 || isempty(Nm), Nm = 31; end
if nargin < 4 || isempty(Nz), Nz = 1001; end

if isscalar(Nm)
  % Compute IFT matrix.
  [IFT,mu] = heateigfun(h0,h0,Nm,Nz);
else
  % IFT and mu as arguments.
  IFT = Nm; mu = Nz;
  Nm = size(IFT,1); Nz = size(IFT,2);
end

% Function that returns temperature at midpoint.
if ~isinf(Nz)
  vmid = @(v) interp1(linspace(0,1,Nz),v,.5);
else
  vmid = @(v) v(.5);
end

if Tcook < .5
  % Find the 1-mode approximation to the cooking time.
  % This works well when the time is comparable to 1/mu(1)^2.
  Teqm = heatsteady(h0,h0,mu.');
  % Eigenfunction evaluated at midpoint.
  IFT1mid = vmid(IFT(1,:));
  tc1 = 1/mu(1)^2*log((Teqm(1)*IFT1mid)/(.5 - Tcook));
else
  % The temperature at z=.5 is less than Tcook, so it's impossible for the
  % food to cook.
  tc = NaN;
  tc1 = NaN;
  return
end

% Set up a function for finding the cooktime using fzero.
f = @(t) vmid(heat(t,0,h0,h0,IFT,mu)) - Tcook;

if Tcook > 1e-2
  tc = fzero(f,[.5*tc1 2*tc1]);
else
  % When Tcook is too small, the 1-mode cooking time differs too much
  % from the actual time to bracket the zero.
  tc = fzero(f,[.5*tc1]);
end
