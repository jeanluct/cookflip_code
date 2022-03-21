function Ts = heatsteady(h0,h1,Nz)
%HEATSTEADY   Steady temperature profile for heat equation.
%   TS = HEATSTEADY(H0,H1) returns the steady-state (equilibrium)
%   temperature profile for the heat equation with unit thermal diffusivity
%   on the interval [0,1], with dimensionless heat transfer coefficients H0
%   at 0 and H1 at 1.  A perfect insulator has H=0, and a perfect conductor
%   has H=inf.  The temperature profile T is a <strong>chebfun</strong> object.
%
%   The boundary conditions at 0 and 1 are
%
%     -DT(0) = H0 (1 - T(0))
%      DT(1) = H1 (0 - T(1))
%
%   where DT denotes the derivative of temperature T.  The temperature
%   difference between the top and the bottom is 1.
%
%   TS = HEATSTEADY(H0,H1,NZ) uses NZ equally-spaced discretization points
%   (default 1001).  Set NZ=Inf for a <strong>chebfun</strong> object as above.
%
%   TM = HEATSTEADY(H0,H1,MU) returns the Fourier components of the
%   steady profile given a vector of eigenvalues MU.
%
%   See also HEATEIGVAL, HEATEIGFUN, CHEBFUN.

% Use "cooking" values if h0, h1 not given.
if nargin < 1 || isempty(h0), h0 = 21.6; end
if nargin < 2 || isempty(h1), h1 = 1.44; end

% Convert h0, h1 to angles.
alpha = acot(h0); beta = acot(h1);

if nargin > 2 && ~isscalar(Nz)
  % Nz contains eigenvalues.
  mu = Nz;
  Cm2 = .5*(1 + (mu./h0).^2) + 1/h0*sin(mu).^2 + ...
        .25*(mu*h0^-2 - mu.^-1).*sin(2*mu);
  Ts = 1./mu./sqrt(Cm2);
  return
end

if nargin < 3 || isempty(Nz), Nz = 1001; end

if isinf(Nz)
  z = chebfun('x',[0 1]);
else
  z = linspace(0,1,Nz).';
end

Ts = cos(alpha)*(sin(beta) + cos(beta)*(1 - z)) / ...
     (sin(alpha + beta) + cos(alpha)*cos(beta));
