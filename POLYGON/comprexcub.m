

function [pts,w,momerr] = comprexcub(deg,X,omega,pos)

%--------------------------------------------------------------------------
% Object:
% Compression of bivariate cubature formulas by Tchakaloff points
% or approximate Fekete points
% useful, for example, in node reduction of algebraic cubature formulas
% see the web page: http://www.math.unipd.it/~marcov/Tchakaloff.html
%--------------------------------------------------------------------------
% Input:
% deg: polynomial exactness degree
% X: 2-column array of point coordinates
% omega: 1-column array of weights
% pos: NNLS for pos=1, QR with column pivoting for pos=0, Fast Lawson-Hanson
%      for pos=2.
%--------------------------------------------------------------------------
% Output:
% pts: 2-column array of extracted points
% w: 1-column array of corresponding weights (positive for pos=1)
% momerr: moment reconstruction error
%--------------------------------------------------------------------------
% Authors:
% Federico Piazzon, Alvise Sommariva and Marco Vianello
% University of Padova, May 2016
%--------------------------------------------------------------------------

% FUNCTION BODY

% troubleshooting
if nargin < 4, pos=3; end

% total-degree Chebyshev-Vandermonde matrix at X
rect=[min(X(:,1)) max(X(:,1)) min(X(:,2)) max(X(:,2))];



V=chebvand(deg,X,rect);
% polynomial basis orthogonalization
[Q,R]=qr(V,0);
% tiny imaginary parts could appear
Q=real(Q);
% possible re-orthogonalization
% [Q,R]=qr(Q,0);

% moments of the orthogonal basis
orthmom=Q'*omega;
% weigths computation
switch pos
    case 1
        % Tchakaloff points (positive weights)
        weights=lsqnonneg(Q',orthmom);
    case 2
        weights = lawsonhanson(Q',orthmom);
    case 3
        weights = LHDM(Q',orthmom);
    otherwise
        % approximate Fekete points (possible negative weights)
        weights=Q'\orthmom;
end
% indexes of nonvanishing weights and compression
ind=find(abs(weights)>0);
pts=X(ind,:);
w=weights(ind);

% moment reconstruction error
% bivariate Chebyshev basis
mom=V'*omega;
momerr=norm(V(ind,:)'*w-mom);
% discrete OP basis
% momerr=norm(Q(ind,:)'*w-orthmom);








%--------------------------------------------------------------------------
% chebvand
%--------------------------------------------------------------------------

function V = chebvand(deg,X,rect)

%--------------------------------------------------------------------------
% Object:
% Vandermonde matrix at X points, relative to Chebyshev basis of degree deg
% on rectangle "rect". The basis is of tensorial type, but of total degree.
%--------------------------------------------------------------------------
% Input:
% deg: polynomial degree.
% X: 2-column array of point coordinates.
% rect: 4-component vector such that the rectangle.
% [rect(1),rect(2)] x [rect(3),rect(4)] contains X.
%--------------------------------------------------------------------------
% Output:
% V: Chebyshev-Vandermonde matrix at X, graded lexic. order.
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

% FUNCTION BODY
% rectangle containing the mesh
if isempty(rect)
    rect=[min(X(:,1)) max(X(:,1)) min(X(:,2)) max(X(:,2))];
end

% couples with length less or equal to deg
% graded lexicographical order
j=(0:1:deg);
[j1,j2]=meshgrid(j);
dim=(deg+1)*(deg+2)/2;
couples=zeros(dim,2);
for s=0:deg
    good=find(j1(:)+j2(:)==s);
    couples(1+s*(s+1)/2:(s+1)*(s+2)/2,:)=[j1(good) j2(good)];
end

% mapping the mesh in the square [-1,1]^2
a=rect(1);b=rect(2);c=rect(3);d=rect(4);
map=[(2*X(:,1)-b-a)/(b-a) (2*X(:,2)-d-c)/(d-c)];

% Chebyshev-Vandermonde matrix on the mesh
T1=chebpolys(deg,map(:,1));
T2=chebpolys(deg,map(:,2));
V=T1(:,couples(:,1)+1).*T2(:,couples(:,2)+1);








%--------------------------------------------------------------------------
% chebpolys
%--------------------------------------------------------------------------

function T=chebpolys(deg,x)

%--------------------------------------------------------------------------
% Object:
% This routine computes the Chebyshev-Vandermonde matrix on the real line
% by recurrence.
%--------------------------------------------------------------------------
% Input:
% deg: maximum polynomial degree
% x: 1-column array of abscissas
%--------------------------------------------------------------------------
% Output:
% T: Chebyshev-Vandermonde matrix at x, T(i,j+1)=T_j(x_i), j=0,...,deg.
%--------------------------------------------------------------------------
% Authors:
% Alvise Sommariva and Marco Vianello
% University of Padova, December 15, 2017
%--------------------------------------------------------------------------

T=zeros(length(x),deg+1);
t0=ones(length(x),1);
T(:,1)=t0;
t1=x;
T(:,2)=t1;

for j=2:deg
    t2=2*x.*t1-t0;
    T(:,j+1)=t2;
    t0=t1;
    t1=t2;
end


