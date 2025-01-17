function [U,jvec,Q,R,dbox,degs] = dORTHVAND(deg,X,u,jvec,C,dbox)

%--------------------------------------------------------------------------
% Object:
% This routine computes a Vandermonde-like matrix for degree "deg" on a 
% d-dimensional point cloud "X" in an total-degree discrete orthogonal 
% polynomial basis w.r.t. the nonnegative weight array "u". 
%--------------------------------------------------------------------------
% Input:
% deg: total degree of the algebraic polynomial basis; 
% X: d-column array of point coordinates;
% * u: 1-column array of nonnegative weights, or nonnegative scalar in
%    case of equal weights;
% * jvec: vector of column indexes that selects a polynomial basis;  
% * C : Chebyshev-Vandermonde matrix on "jvec" basis
% * dbox: variable that defines a hyperrectangle with sides parallel to the
%    axis, containing the domain (or "X" in the discrete case). 
%    If "dbox" is not provided, it is the smaller "hyperrectangle", with 
%    sides parallel to the cartesian axes, containing the pointset "X".
%    It is a matrix with dimension "2 x d", where "d" is the dimension of
%    the space in which it is embedded the domain. 
%    For instance, for a 2-sphere, it is "d=3", for a 2 dimensional polygon 
%    it is "d=2".
%    As example, the set "[-1,1] x [0,1]" is described as "[-1 0; 1 1]".
% * dimpoly: dimension of polynomial space (if known in advance, useful
%   for instance on the sphere or its portions). If not available, set 
%   " dimpoly='' " or do not declare the variable.
% Note: the variables with an asterisk "*" are not mandatory and can be 
% also set as empty matrix.
%--------------------------------------------------------------------------
% Output:
% U: Vandermonde-like matrix in a "u"-orthogonal polynomial basis on "X"; 
% jvec: vector of column indexes, selects a polynomial basis;
% Q: orthogonal factor in the QR decomposition
%                  diag(sqrt(u))*C(:,jvec)=Q*R 
%    where "C=dCHEBVAND(n,X,dbox)" (if not assigned in inputs);
% R: triangular matrix of the QR decomposition
% dbox: it defines a hyperrectangle with sides parallel to the axis, 
%    containing the domain (or the pointset if "dbox" was not assigned as
%    input variable.
%--------------------------------------------------------------------------
% Dates:
%
% Written on 26/07/2020 by M. Dessole, F. Marcuzzi, M. Vianello.
% Modified on: 21/12/2022: A.Sommariva.
%--------------------------------------------------------------------------



% ........................... Function body ...............................



% ...... troubleshooting ......
if nargin < 7, dimpoly=[]; end
if nargin < 6, dbox=[]; end
if nargin < 5, C=[]; end
if nargin < 4, jvec=[]; end
if nargin < 3, u=[]; end

if isempty(C), [C,degs,dbox0] = chebvand(deg,X,dbox); end
if isempty(dbox), dbox=dbox0; end
if isempty(jvec), N=rank(C); else, N=length(jvec); end
if isempty(u), u=1; end

% ......................... main code below ...............................

% ...... computing Vandermonde matrix "U" ......

% scaling the matrix rows by the sqrt of the weights 
B = zeros(size(C));
for k=1:length(C(1,:))
    B(:,k)=C(:,k).*sqrt(u);
end

% polynomial basis orthogonalization
if N < length(C(1,:)) % low rank (e.g. the domain is an algebraic surface)
    if isempty(jvec)
        [Q0,R0,pm]=qr(B,0); 
        jvec=pm(1:N);
        R=R0(1:N,1:N);
        Q=Q0(:,1:N);
    else
        [Q,R]=qr(B(:,jvec),0);
    end
else % full rank
    [Q,R]=qr(B,0);
    jvec=(1:N);
end

% evaluation of the orthogonal matrix
U=C(:,jvec)/R;




