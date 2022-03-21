function [varargout] = mincooktime(Nflips,Tcook,h0,h1,Nm,Nz)
%MINCOOKTIME   Minimize cooking time for given number of flips of the food.
%   [TC,TF] = MINCOOKTIME(NFLIP,COOKTEMP,...) returns the total cooking time
%   TC and the optimal flip times TF for exactly NFLIP flips of the food.
%   The value of COOKTEMP and following parameters are passed to <strong>cooktime</strong>.
%
%   See also COOKTIME, HEATEIGFUN, HEATEIGVAL. HEATSTEADY.

if nargin < 1 || isempty(Nflips), Nflips = 1; end
if nargin < 2 || isempty(Tcook), Tcook = .257; end
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

% The asymptotic total cooktime for infinite flips.
% Eventually use a theoretical value as an initial guess.
if isinf(h0) && isinf(h1)
  % For the symmetric conducting case, the optimal time is essentially
  % constant as a function of the number of flips.
  tcookinf = .0976;
elseif h0 == 2 && h1 == 2
  tcookinf = .281725;
else
  % This is an extrapolated value for the default h0,h1.
  tcookinf = .0754;
end

% Start all intervals equal.
tstart = repmat(tcookinf/(Nflips+1),[1 Nflips]);

f = @(tflips) cooktime(tflips,Tcook,h0,h1,IFT,mu);

if nargout > 0
  opts = optimset('Display','off');
else
  opts = optimset('Display','iter');
end
%opts = optimset(opts,'TolX',1e-8);

[tflipsopt,ttotalopt] = fminsearch(f,tstart,opts);

if nargout > 0
  varargout{1} = ttotalopt;
  if nargout > 1
    varargout{2} = tflipsopt;
  end
else
  % Call f again to make plot.
  f(tflipsopt)
end
