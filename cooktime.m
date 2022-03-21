function [varargout] = cooktime(tflip,Tcook,h0,h1,Nm,Nz)
%COOKTIME   Cooking time on a hot plate for several flips of the food.
%   TC = COOKTIME(TFLIP,COOKTEMP) returns the total cooking time TC for a
%   one-dimensional slab of food.  The food is considered cooked when every
%   point has experienced a temperature greater than COOKTEMP.  TFLIP is a
%   vector containing the duration of each "flip" of the food on the hot
%   plate.  Set TFLIP=inf to never flip and get the "cookthrough" time,
%   if it is defined.
%
%   TC = COOKTIME(TFLIP,COOKTEMP,H0,H1) uses dimensionless heat transfer
%   coefficients H0 (hot plate) and H1 (free top surface).
%
%   If the food never cooks in the allotted time, TC is returned as NaN.  If
%   the food cooks before the final flip, TC is also returned as NaN.
%
%   TC = COOKTIME(TFLIP,COOKTEMP,H0,H1,NM,NZ) uses NM eigenfunctions with NZ
%   spatial gridpoints.
%
%   TC = COOKTIME(TFLIP,COOKTEMP,H0,H1,IFT,MU) uses the matrix IFT and
%   vector MU returned by <strong>heateigfun</strong> rather than computing
%   them from scratch.
%
%   [TC,t,COOKFRAC,z,TMAX] = COOKTIME(...) returns the time series for the
%   cooked fraction COOKFRAC as a function of time t, as well as the profile
%   of maxmimum temperatures achieved TMAX as a function of z.
%
%   See also HEATEIGFUN, HEATSTEADY, MINCOOKTIME, TCOOKTHRU.

% Use optimal solution as defaults.
if nargin < 1 || isempty(tflip), tflip = 0.045625; end
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

% Plot steps (time data) when there are no output arguments.
% This is so mincooktime can call cooktime many times, rapidly.
plotsteps = (nargout < 1);
% But if we demand the time data, then make it finer.
% (Presumably then we're not being called by mincooktime.)
finedata = (plotsteps || nargout > 1);

if exist('plotprops') == 2
  plotprops
else
  % Font and graphics properties for pretty plots.
  fonttype = 'Times';
  fsize = 16; fcsize = 12; lw = 2; ms = 15;
  txtattrib = {'FontName',fonttype,'FontSize',fsize,'FontWeight','normal'};
  txtattribtex = {txtattrib{:},'Interpreter','Latex'};
end

m = 1:Nm;

if ~isinf(Nz)
  z = linspace(0,1,Nz);
  intcols = @(f) trapz(z,f,2);     % integral of columns
else
  % Unfortunately using chebfun turns out to be prohibitively slow.
  z = chebfun('x',[0 1]);
  intcols = @(f) integral(f);      % integral of columns
end

% Equilibrium profile and its Fourier coefficients.
Teq = heatsteady(h0,h1,Nz); Teqm = heatsteady(h0,h1,mu.');

if Tcook > Teq(1)
  error('Tcook cannot be greater than T(z=0)=%g.',Teq(1));
end

if isinf(tflip(1))
  % If we're not flipping, find the cookthrough time.
  % The time required to cook without flipping.
  if isnan(tcookthru(Tcook,h0,h1,IFT,mu))
    error('Flip at least once to cook since Tcook > T(z=1)=%g.',Teq(end))
  end
end

% "Flipping" operation.
% The "flipping" matrix: project flipped modes onto unflipped modes.
f = flipop(IFT);

% The flipped equilibrium profile.
B = Teqm - f*Teqm;

% Initial temperature deviation = IFT of Teq.
thetam0 = -Teqm;
thetam = thetam0;

% The number of integration steps to take per flipping interval.
% Use relaxation time to determine baseline step size dt0.
% TODO: select these tolerances on command-line?
trelax = 1/mu(1)^2; tlast = 100*trelax;
if plotsteps || finedata
  % If we're plotting, use a uniform step for good resolution, but not
  % particularly accurate cooking time.
  dt0 = trelax/300;
  dtlast = dt0;
else
  % If we're not plotting, use a coarse step before the last flip, then
  % slow down to get an accurate estimate.
  dt0 = trelax/10;
  dtlast = trelax/1000;
end
% Compute number of steps based on dt0.
if ~isinf(tflip(1))
  Nsteps = [ceil(tflip/dt0) ceil(tlast/dtlast)];
  % Now correct dt in each interval to ensure we land at exact boundaries.
  dt = [tflip tlast]./Nsteps;
else
  Nsteps = ceil(tlast/dtlast);
  % Now correct dt in each interval to ensure we land at exact boundaries.
  dt = tlast./Nsteps;
end

% Points are "cooked" if they achieve the temperature Tcook at any time.
cookedfrac = [0];
zzl = [NaN NaN];  % list of zeros of Tmax; should be at most two
if ~isinf(Nz)
  Tmax = zeros(Nz,1);
else
  Tmax = 0*z;
end

% Midpoint temperature.
Tmid = [0];

if plotsteps, figure(1); end

if ~isinf(tflip(1))
  jmax = length(tflip)+1;
else
  jmax = 1;
end

for j = 1:jmax
  % Cook for a time interval tflip(j).
  for i = 1:Nsteps(j)
    t = sum(tflip(1:j-1)) + i*dt(j);

    % Evolve temperature for a time dt.
    thetam = thetam.*exp(-mu.^2*dt(j))';
    T = Teq + IFT'*thetam;

    % Find the maximum temperature achieved at each point.
    % This operation really kills chebfun.
    Tmax = max(T,Tmax);
    [cf,zz] = cookedfraction(z,Tmax - Tcook);
    cookedfrac = [cookedfrac cf];
    zzl = [zzl ; zz];

    % Temperature at midpoint.
    Tmid = [Tmid T(ceil(Nz/2))];

    if plotsteps
      plot(T,z,'b-','LineWidth',lw)
      hold on
      if false
        %plot(flipud(T),z,'b--','LineWidth',lw)
        % Plot antisymmetric part of profile, centered at .5 for convenience.
        plot(.5 - .5*(T - flipud(T)),z,'m-','LineWidth',lw)
        %plot(.5*(T + flipud(T)),z,'m-','LineWidth',lw)
      end
      if j-1 ~= 1, s = 's'; else s = []; end
      title(sprintf(['%g flip' s],j-1),txtattribtex{:})
      text(.7,.6,sprintf('t=%.3f',t))
      text(.7,.5,sprintf('%3.0f%% cooked',100*cookedfrac(end)))
      hold on
      if false
        % Plot the equilibrium profile.
        plot(Teq,z,'k--','LineWidth',.5*lw)
      end
      plot(Tcook*[1 1],[0 1],'r--','LineWidth',.5*lw)
      if ~isinf(Nz)
        Tover = T.*(Tmax >= Tcook); Tover(Tover == 0) = NaN;
      else
        Tover = (Tmax >= Tcook);
        Tover = T.*Tover + (1 - Tover)*(-100);
      end
      plot(Tover,z,'r-','LineWidth',2*lw)
      hold off
      xlabel('$T$',txtattribtex{:})
      ylabel('$z$',txtattribtex{:})
      set(gca,txtattrib{:})
      axis([0 1 0 1])
      drawnow
    end

    if cookedfrac(end) == 1
      if ~nargout
        [~,ilastcooked] = min(Tmax);
        zlastcooked = z(ilastcooked);
        fprintf('Cooked!  Last point cooked was z=%.5f\n',zlastcooked);
      end
      break
    end
    if plotsteps && i == Nsteps(j), pause; end
  end

  % break only breaks out of innermost loop.  Must break again if cooked.
  if cookedfrac(end) == 1, break; end

  % Flip!
  thetam = f*thetam - B;
  Tmax = flipud(Tmax);
end

% If uncooked, return NaN.
if cookedfrac(end) < 1
  warning('cooktime:uncooked','Uncooked after %d flips.',length(tflip))
  ttotal = NaN;
  if nargout > 0, varargout{1} = ttotal; end
  return
end

% If cooked before final interval, return NaN as well, since this
% indicates the subsequent flips were not needed.
if j < length(tflip)+1 && ~isinf(tflip(1))
  warning('cooktime:cookedbeforefinal','Cooked before final flip.')
  ttotal = NaN;
  if nargout > 0, varargout{1} = ttotal; end
  return
end

% Adjust Nsteps and tflip in final interval to reflect the steps actually taken.
%
% TODO: Note this assumes the final interval was needed.  The food may cook
% completly on an earlier interval.  Handle this better.
Nsteps(j) = i;
if ~isinf(tflip(1))
  tflip = [tflip dt(j)*Nsteps(j)];
else
  tflip = dt(j)*Nsteps(j);
end
t = cumsum(tflip);
ttotal = t(end);

trange = [0];
for j0 = 1:j
  trange = [trange trange(end)+dt(j0)*(1:Nsteps(j0))];
end

if nargout > 0
  varargout{1} = ttotal;
  if nargout > 1
    varargout{2} = trange;
  end
  if nargout > 2
    varargout{3} = cookedfrac;
  end
  if nargout > 3
    varargout{4} = z;
  end
  if nargout > 4
    varargout{5} = Tmax;
  end
  if nargout > 5
    % Undocumented: return zeros of Tmax.
    % zzl is a vector of two columns containing <= 2 zeros.
    % If there is only one zero, the second column is NaN.
    % If there are no zeros, both columns are NaN.
    % The first and last entry are both [NaN NaN], since at first no points
    % are cooked, and at the end all points are cooked.
    varargout{6} = zzl;
  end
else
  fmt0 = ' = %.5f'; fmt1 = '  (%.2f%%)'; fmt = [fmt0 fmt1 '\n'];
  if length(tflip) > 1
    for j0 = 1:j-1
      fprintf(['  dt%-2d' fmt],j0,tflip(j0),100*tflip(j0)/ttotal)
    end
    fprintf(['dtlast' fmt],tflip(j),100*tflip(j)/ttotal)
  end
  fprintf(['ttotal' fmt],ttotal,100)
  if length(tflip) > 2
    fprintf([' dtavg' fmt0],mean(tflip(1:j-1)))
    fprintf([' +/- %.5f (not including dtlast)\n'],std(tflip(1:j-1)))
  end

  figure(2)
  plot(trange,cookedfrac,'k-','LineWidth',lw)
  xlabel('$t$',txtattribtex{:})
  ylabel('cooked fraction',txtattribtex{:})
  if length(tflip)-1 ~= 1, s = 's'; else s = []; end
  title(sprintf(['%d flip' s],length(tflip)-1),txtattribtex{:})
  set(gca,txtattrib{:})
  axis tight
  hold on
  tlast = 0; Nstepsc = cumsum(Nsteps);
  xl = xlim; xr = xl(2)-xl(1); xlabpos = .05*xr;
  yl = ylim; yr = yl(2)-yl(1); ylabpos = .05*yr;
  for j0 = 1:j-1
    % The cooked fraction at the end of t(j0).
    % The '+1' is due to 0 being recorded at the start.
    cookedfrac1 = cookedfrac(Nstepsc(j0)+1);
    % Plot dashed lines indicating the cooked and time fraction at the end of
    % tflip(j).
    plot([0 t(j0)],cookedfrac1*[1 1],'k--','LineWidth',lw)
    plot(t(j0)*[1 1],cookedfrac1*[0 1],'k--','LineWidth',lw)
    if length(tflip) < 14
      %text(t(j0)/4,cookedfrac1+.05,sprintf('%3.0f%% cooked',100*cookedfrac1))
      text(xlabpos,cookedfrac1+.03*yr,sprintf('%3.0f%%',100*cookedfrac1))
      %text(1.02*t(j0),cookedfrac1/3,sprintf('%3.0f%% time',100*t(j0)/ttotal))
      text(t(j0)+.007*xr,ylabpos,sprintf('%3.0f%%',100*t(j0)/ttotal))
    end
  end
  hold off

  figure(3)
  plot(trange,Tmid,'k-','LineWidth',lw)
  xlabel('$t$',txtattribtex{:})
  ylabel('midpoint temperature $T(1/2)$',txtattribtex{:})
  if length(tflip)-1 ~= 1, s = 's'; else s = []; end
  title(sprintf(['%d flip' s],length(tflip)-1),txtattribtex{:})
  set(gca,txtattrib{:})
  axis tight
  hold off

  figure(4)
  plot(Tmax,z,'b-','LineWidth',lw)
  if j-1 ~= 1, s = 's'; else s = []; end
  title(sprintf(['max temperature after %g flip' s],j-1),txtattribtex{:})
  hold on
  plot(Tcook*[1 1],[0 1],'r--','LineWidth',.5*lw)
  axis([0 1 0 1])
  hold off
  xlabel('$T_{\mathrm{max}}$',txtattribtex{:})
  ylabel('$z$',txtattribtex{:})
  set(gca,txtattrib{:})
end

%============================================================================
function [cf,zz] = cookedfraction(z,dT)
%COOKEDFRACTION   Compute the fraction of food cooked.

zz = findzeros(z,dT);  % find zeros of Tmax - Tcook

Nz = max(size(dT));
if ~isinf(Nz)
  dT0 = dT(1);
else
  dT0 = dT(0);
end

if isempty(zz)
  % If there are no zeros, we're either all cooked or not cooked at all.
  if dT0 >= 0, cf = 1; else cf = 0; end
  zz = [NaN NaN];  % no zeros, so fill with NaN
elseif length(zz) == 1
  % If there's one zero (say, first phase), then the bottom or top half
  % is cooked.
  if dT0 >= 0
    cf = zz;
    zz = [zz NaN];  % bottom is cooked (should always be the case)
  else
    cf = 1-zz;
    zz = [NaN zz];  % top is cooked (weird?)
    warning('Only top is cooked?')
  end
elseif length(zz) == 2
  % If there's two zeros, there are three cooked intervals.
  if dT0 >= 0
    cf = 1 - (zz(2) - zz(1));
  else
    cf = zz(2) - zz(1);
  end
else
  % If there are more zeros, something is weird.
  warning('More than two zeros?')
end

% Trust but verify.
if ~isinf(Nz)
  cf_coarse = @(dT) nnz(dT >= 0)/Nz;
  if abs(cf - cf_coarse(dT)) > 2/Nz
    warning('Inconsistency in cooked fraction (%g vs %g).',cf,cf_coarse(dT))
  end
end

%============================================================================
function z = findzeros(x,f)
%FINDZEROS   Find all zeros of the function f.

if ~isinf(max(size(f)))
  z = [];

  sc = sign(f(2:end)) - sign(f(1:end-1));

  for i = 1:length(sc)
    if abs(sc(i)) == 1
      z = [z x(i)];
    elseif abs(sc(i)) == 2
      ii = [i i+1];
      z = [z interp1(f(ii),x(ii),0)];
    end
  end

  if ~nargout
    plot(x,f)
    hold on
    plot(z,zeros(size(z)),'r.')
    hold off
  end

else
  z = roots(f).';
end
