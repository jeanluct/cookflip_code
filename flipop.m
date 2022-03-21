function f = flipop(h0,h1,Nm,Nz)
%FLIPOP   The "flipping" operator acting on eigenfunctions.
%   F = FLIPOP(H0,H1,NM,NZ) returns the "flipping" operator, which is
%   defined as
%
%     F(m,n) = < V_m(1-z) , V_n(z) >
%
%   where < , > is the inner product, with |V_n|^2 = 1.  See
%   <strong>heateigfun</strong> for a description of H0,H1,NM,NZ.
%
%   F = FLIPOP(IFT) uses the matrix IFT returned by <strong>heateigfun</strong>
%   rather than computing it from scratch.
%
%   See also HEATEIGFUN, FLIPHEATFIX.

%
% This file is part of cookflip_code
%
% Copyright (c) 2022 Jean-Luc Thiffeault <jeanluc@math.wisc.edu>
%
% See the file LICENSE for copying permission.
%

% Use "cooking" values if h0, h1 not given.
if nargin < 1 || isempty(h0), h0 = 21.6; end
if nargin < 2 || isempty(h1), h1 = 1.44; end
if nargin < 3 || isempty(Nm), Nm = 20; end
if nargin < 4 || isempty(Nz), Nz = 1001; end

if isscalar(h0)
  % Compute IFT matrix.
  IFT = heateigfun(h0,h1,Nm,Nz);
else
  % First arg is IFT matrix.
  IFT = h0;
  Nm = size(IFT,1); Nz = size(IFT,2);
end

if ~isinf(Nz)
  z = linspace(0,1,Nz);
  intcols = @(f) trapz(z,f,2);     % integral of columns
else
  intcols = @(f) integral(f);      % integral of columns
end

% The "flipping" matrix: project flipped modes onto unflipped modes.
f = zeros(Nm,Nm);
for i = 1:Nm
  for j = 1:Nm
    f(i,j) = intcols(IFT(i,:).*fliplr(IFT(j,:)));
  end
end
