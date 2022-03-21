function [tct,tct1] = tcookthru(Tcook,h0,h1,Nm,Nz)
%TCOOKTHRU   The time to cook food through without flipping.
%   TCT = TCOOKTHRU(COOKTEMP,H0,H1,NM,NZ) returns the time TCT needed to cook a
%   slab of food of unit thickness, without flipping.  TCT=NaN if this is
%   not possible (i.e., the top temperature is below COOKTEMP).  See
%   <strong>heateigfun</strong> for a description of H0,H1,NM,NZ.
%
%   TCT = TCOOKTHRU(COOKTEMP,H0,H1,IFT,MU) uses the matrix IFT and vector MU
%   returned by <strong>heateigfun</strong> rather than computing them from scratch.
%
%   [TCT,TCT1] = TCOOKTHRU(...) also returns TCT1, the 1-mode approximation
%   to the cookthrough time.
%
%   See also HEATEIGFUN, TEMP.

%
% This file is part of cookflip_code
%
% Copyright (c) 2022 Jean-Luc Thiffeault <jeanluc@math.wisc.edu>
%
% See the file LICENSE for copying permission.
%

if nargin < 1 || isempty(Tcook), Tcook = .257; end
% Use "cooking" values if h0, h1 not given.
if nargin < 2 || isempty(h0), h0 = 21.6; end
if nargin < 3 || isempty(h1), h1 = 1.44; end
if nargin < 4 || isempty(Nm), Nm = 31; end
if nargin < 5 || isempty(Nz), Nz = 1001; end

if isscalar(Nm)
  % Compute IFT matrix.
  [IFT,mu] = heateigfun(h0,h1,Nm,Nz);
else
  % IFT and mu as arguments.
  IFT = Nm; mu = Nz;
  Nm = size(IFT,1); Nz = size(IFT,2);
end

Teq = heatsteady(h0,h1,Nz); Teqm = heatsteady(h0,h1,mu.');

if Teq(end) > Tcook
  % Find the 1-mode approximation to the cookthrough time.
  % This works well when the time is comparable to 1/mu(1)^2.
  tct1 = 1/mu(1)^2*log((Teqm(1)*IFT(1,end))/(Teq(end) - Tcook));
else
  % The temperature at z=1 is less than Tcook, so it's impossible for the
  % food to cook through.
  tct = NaN;
  tct1 = NaN;
  return
end

% Set up a function for finding the cooktime using fzero.
vend = @(v) v(end);
f = @(t) vend(heat(t,0,h0,h1,IFT,mu)) - Tcook;

if Tcook > 1e-2
  tct = fzero(f,[.5*tct1 2*tct1]);
else
  % When Tcook is too small, the 1-mode cooking time differs too much
  % from the actual time to bracket the zero.
  tct = fzero(f,[.5*tct1]);
end
