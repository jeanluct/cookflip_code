function [varargout] = heateigfun(h0,h1,mu,Nz)
%HEATEIGFUN   Eigenfunctions of the heat operator.
%   [IFT,MU] = HEATEIGFUN(H0,H1,N,NZ) returns the first N eigenfunctions of
%   the heat operator with unit thermal diffusivity on the interval [0,1],
%   with dimensionless heat transfer coefficients H0 at 0 and H1 at 1.  A
%   perfect insulator has H=0, and a perfect conductor has H=inf.  The
%   corresponding (square root of) eigenvalues are returned in the array MU,
%   and are ordered by increasing size.
%
%   The L2-normalized eigenfunctions IFT are returned as the rows of an N x
%   NZ array, with NZ equally-spaced discretization points.  If NZ=Inf, a
%   <strong>chebfun</strong> quasimatrix object of size N x Inf is returned instead.
%   Set NZ=[] to use default NZ=1001.
%
%   IFT can be used to take an "inverse Fourier transform" of a
%   column-vector TM of length N by taking the product TZ = TM'*IFT.
%
%   IFT can also be used to take the Fourier transform of a function TZ
%   represented by a vector of size NZ by taking TM = IFT*TZ/(NZ-1).  For
%   NZ=Inf, use TM = IFT*TZ since <strong>chebfun</strong> then carries out
%   the integral accurately.
%
%   IFT = HEATEIGFUN(H0,H1,MU,NZ) returns the eigenfunctions given
%   precomputed eigenvalues MU, as returned by <strong>heateigval</strong>.
%
%   References
%
%   T. A. Driscoll, N. Hale, and L. N. Trefethen, editors, Chebfun Guide,
%   Pafnuty Publications, Oxford, 2014.  http://www.chebfun.org/docs/guide/
%
%   "Quasimatrices and Least-Squares,"
%   http://www.chebfun.org/docs/guide/guide06.html
%
%   See also HEATEIGVAL, HEATSTEADY, CHEBFUN.

% Use "cooking" values if h0, h1 not given.
if nargin < 1 || isempty(h0), h0 = 21.6; end
if nargin < 2 || isempty(h1), h1 = 1.44; end

if nargin < 3 || isempty(mu)
  % If mu not given, compute 20 modes.
  mu = heateigval(h0,h1,20);
elseif isscalar(mu)
  % If mu is a scalar, use mu modes.
  mu = heateigval(h0,h1,mu);
end

% Use chebfun by default (Nz = inf).
if nargin < 4 || isempty(Nz), Nz = 1001; end

Nm = length(mu);
m = 1:Nm;

if ~isinf(Nz)
  z = linspace(0,1,Nz);
else
  z = chebfun('x',[0 1]);
end

% Make plots when there is no output argument.
ploteigf = (nargout < 1);

% Convert h0, h1, to angles.
alpha = acot(h0); beta = acot(h1);

%
% Eigenfunctions:
%
% Can be used as operator for (slow) inverse Fourier transform:
%
% Given a row-vector of mode coefficients Tm, transform to z-space by
%
% T = Tm*IFT
%
if ~isinf(Nz)
  IFT = cos(alpha)*sin(kron(mu.',z)) + sin(alpha)*diag(mu)*cos(kron(mu.',z));
  % Normalize eigenfunctions.
  intcols = @(f) trapz(z,f,2);     % integral of columns
  %intcols = @(f) sum(f,2)/(Nz-1);  % less accurate, but then we can compute FT
                                    % by multiplication.
  C2 = intcols(IFT.^2);
  IFT = diag(sqrt(C2))\IFT;
else
  % For Nz=Inf, use chebfun.
  IFT = [];
  for i = m
    IFT = [IFT cos(alpha)*sin(mu(i)*z) + sin(alpha)*mu(i)*cos(mu(i)*z)];
  end
  IFT = IFT * diag(1./sqrt(integral(IFT.^2))); % normalize
  IFT = IFT.';
  intcols = @(f) integral(f);      % integral of columns
end

if nargout > 0
  varargout{1} = IFT;
  if nargout > 1
    varargout{2} = mu;
  end
end

if ploteigf
  % Plot the first 5 eigenfunctions.
  for j = 1:min(5,length(m))
    plot(z,IFT(j,:))
    hold on
  end
  hold off

  ortho = zeros(Nm,Nm);
  for i = m
    for j = m
      ortho(i,j) = intcols(IFT(i,:).*IFT(j,:));
    end
  end
  orthoerr = max(max(abs(ortho - eye(Nm))));
  fprintf('Error in orthonormality = %.2g\n',orthoerr)
  if orthoerr > 1e-3
    warning('Modes don''t appear orthonormal (error = %.2e).',orthoerr)
  end
end
