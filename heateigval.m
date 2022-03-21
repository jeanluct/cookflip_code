function [varargout] = heateigval(h0,h1,Nm)
%HEATEIGVAL   Eigenvalues (spectrum) of the heat operator.
%   MU = HEATEIGVAL(H0,H1,N) returns the square root MU of the N smallest
%   eigenvalues of the heat operator with unit thermal diffusivity on the
%   interval [0,1], with dimensionless heat transfer coefficients H0 at 0
%   and H1 at 1.  A perfect insulator has H=0, and a perfect conductor has
%   H=inf. The eigenvalues are ordered by increasing size.
%
%   Eigenfunctions of the heat equation decay in time at the rate
%   exp(-MU.^2*t).
%
%   See also HEATEIGFUN, HEATSTEADY.

% Use "cooking" values if h0, h1 not given.
if nargin < 1 || isempty(h0), h0 = 21.6; end
if nargin < 2 || isempty(h1), h1 = 1.44; end
if nargin < 3 || isempty(Nm), Nm = 20; end

%
% TODO:
%
% - Cases with negative eigenvalues?
%

if h0 < 0 || h1 < 0
  error('Only positive h0 and h1 supported for now.')
end

% Make plots when there is no output argument.
ploteigv = (nargout < 1);

% Convert h0, h1 to angles.
alpha = acot(h0); beta = acot(h1);
sab = sin(alpha + beta);
sasb = sin(alpha)*sin(beta); cacb = cos(alpha)*cos(beta);

% The eigenvalue equation.
% Divide by 1 + sasb*mu.^2 to avoid large values at large mu.
transc = @(mu) (sab*mu.*cos(mu) + (cacb - sasb*mu.^2).*sin(mu)) ./ ...
         (1 + sasb*mu.^2);

opts = optimset('Display','off','TolX',eps);

% There can be two roots in the interval where muc falls.
% As far as I can see this is the only time an interval [2n-1,2n+1]pi/2
% can have more than one root in it.
muc = sqrt(cacb/sasb);
% The intervals that bracket roots.
ivals = unique([0 muc (2*(1:Nm+1)-1)*pi/2]);
mu = zeros(1,Nm);
neigs = 0; nint = 1;
while neigs < Nm
  ival = [ivals(nint) ivals(nint+1)];
  if ploteigv
    figure(1)
    mu0 = linspace(ival(1),ival(2),100);
    plot(mu0,transc(mu0),'k-')
    hold on
  end
  if alpha == 0 && beta == pi/2 || alpha == pi/2 && beta == 0
    % For the zero temp/insulating case, the eigenvalues are odd multiples
    % of pi/2, which are the boundaries of our search intervals.
    % Just solve this problem separately.
    neigs = neigs + 1;
    mu(neigs) = (2*nint-1)*pi/2;
    if ploteigv
      plot(mu(neigs),0,'r.','MarkerSize',15)
    end
  elseif prod(transc(ival)) < 0
    [mu0,fval,exitflag] = fzero(transc,ival,opts);
    if exitflag ~= 1
      error('fzero returned exitflag=%g.\n',exitflag)
    end
    neigs = neigs + 1;
    mu(neigs) = mu0;
    if ploteigv
      plot(mu(neigs),0,'r.','MarkerSize',15)
    end
  end
  nint = nint + 1;
end

if length(mu) ~= Nm
  error('Something went wrong... found only %d < %d modes.',length(mu),Nm)
end

if nargout > 0
  varargout{1} = mu;
end

if ploteigv
  figure(1), hold off;
  fprintf('Max error in eigenval eqn: %.2e\n',max(abs(transc(mu))))
end
