function u = flipheatbl(Z,h0,h1,dt)
%FLIPHEATBL   Boundary layer for the flip-heat operator for rapid flip time.
%   U = FLIPHEATBL(Z,H0,H1) returns the rescaled boundary layer profile U(Z)
%   near Z=0 for the fixed point of the flip-heat operator, in the limit as
%   the flipping interval T goes to zero.  The heat transfer coefficients
%   are H0 and H1.
%
%   FLIPHEATBL(Z) uses defaults H0 = H1 = inf.
%
%   See also FLIPOP, FLIPHEATOP, FLIPHEATFIX.

%
% This file is part of cookflip_code
%
% Copyright (c) 2022 Jean-Luc Thiffeault <jeanluc@math.wisc.edu>
%
% See the file LICENSE for copying permission.
%

if nargin < 1 || isempty(Z), Z = linspace(0,8,100); end
if nargin < 2 || isempty(h0), h0 = inf; end
if nargin < 3 || isempty(h1), h1 = h0; end

if isinf(h0) & isinf(h1)
  fint = @(k,Z) -2/pi*(1 + exp(k.^2)).^-1.*sin(k*Z)./k;

  u = zeros(size(Z));

  for i = 1:length(u)
    fintk = @(k) fint(k,Z(i));
    u(i) = integral(fintk,0,inf);
  end

  u = u + 1;
else
  % Kludge: since we don't have a good integral formula for cases
  % other than h0=h1=inf, compute the boundary layer solution by
  % rescaling flipheatfix.
  % Need to make sure to use dt small enough and enough modes.
  % Amplitude depends on dt.
  if nargin < 4 || isempty(dt), dt = .001; end
  Npts = ceil(length(Z)/sqrt(dt));
  Unum = flipheatfix(dt,h0,h1,[],Npts);
  z = linspace(0,1,Npts);
  u = interp1(z,Unum,sqrt(dt)*Z);
end
