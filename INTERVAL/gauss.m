% GAUSS Gauss quadrature rule.
%
%    Given a weight function w encoded by the nx2 array ab of the 
%    first n recurrence coefficients for the associated orthogonal
%    polynomials, the first column of ab containing the n alpha-
%    coefficients and the second column the n beta-coefficients, 
%    the call xw=GAUSS(n,ab) generates the nodes and weights xw of
%    the n-point Gauss quadrature rule for the weight function w.
%    The nodes, in increasing order, are stored in the first 
%    column, the n corresponding weights in the second column, of
%    the nx2 array xw.
%  
%    See the reference Gaussian Quadrature and the Eigenvalue Problem
%    Theorem 11.
function xw=gauss(N,ab)
N0=size(ab,1); if N0<N, error('input array ab too short'), end
J=zeros(N);
for n=1:N, J(n,n)=ab(n,1); end % The diagonal elements of J are determined
for n=2:N
  J(n,n-1)=sqrt(ab(n,2)); % The lower diagonal elements of J are determined
  J(n-1,n)=J(n,n-1); % The upper diagonal elements of J are determined
end

% J is a triangle symmertric diagonal matrix
[V,D]=eig(J);

% Sort every coloumn in ascending
[D,I]=sort(diag(D)); 

% Resort V by the new order of I
V=V(:,I);

xw=[D ab(1,2)*V(1,:)'.^2];
