function [HF,HFH] = flipheatop(t,h0,h1,Nm,Nz)
%FLIPHEATFOP   The flip-heat operator.
%   HF = FLIPHEATOP(t,H0,H1,NM,NZ) returns the the flip-and-heat operator
%   acting on eigenfunctions of the heat operator for a time t.  HF is a
%   NMxNm matrix.  See <strong>heateigfun</strong> for a description of H0,H1,NM,NZ.
%
%   TF = FLIPHEATOP(t,F,MU) uses the flip operator F and eigenvales MU
%   returned by <strong>fliop</strong> and <strong>heateigval</strong> rather than computing them from scratch.
%
%   [HF,HFH] = FLIPHEATOP(...) also returns the symmetric heat-flip-heat operator.
%
%   See also FLIPOP, FLIPHEATFIX, HEATEIGVAL, HEATEIGFUN.

if nargin < 1 || isempty(t), t = .02; end
% Use "cooking" values if h0, h1 not given.
if nargin < 2 || isempty(h0), h0 = 21.6; end
if nargin < 3 || isempty(h1), h1 = 1.44; end
if nargin < 4 || isempty(Nm), Nm = 31; end
if nargin < 5 || isempty(Nz), Nz = 1001; end

if isscalar(h0)
  % Compute IFT matrix.
  [IFT,mu] = heateigfun(h0,h1,Nm,Nz);
  % The "flipping" matrix: project flipped modes onto unflipped modes.
  F = flipop(IFT);
else
  % F and mu as arguments.
  F = h0; mu = h1;
end

% The flip-heat operator.
HF = diag(exp(-mu.^2*t))*F;

if nargout > 1
  % The symmetric heat-flip-heat operator.
  H = diag(exp(-mu.^2*t/2));
  HFH = H*F*H;
end
