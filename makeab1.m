function [A,b] = makeab1(K,N,Dx,phi0L,phi0R)
%[A,b] = makeab1(K,N,Dx,phi0L,phi0R)
%creates the m by m matrix A and m by 1 vector b of forward solver
%where:
%  K is m by 1 vector of conductances
%  N is m by 1 vector of recharge rates
%  Dx is size of blocks (same for all)
%  phi0L and phi0R are left and right boundary head values.  
%Note that A is sparse.  

%pkk  May 5, 2009

m = length(K);
A = sparse(zeros(m,m)); b = zeros(m,1);
%Mass balance for first block:
A(1,1) = -2*(K(1)+K(1)*K(2)/(K(1)+K(2)))/Dx;
A(1,2) = 2*K(1)*K(2)/(K(1)+K(2))/Dx;
b(1) = -N(1)*Dx-2*phi0L*K(1)/Dx;
%Mass balance for last block:
A(m,m) = -2*(K(m)+K(m)*K(m-1)/(K(m)+K(m-1)))/Dx;
A(m,m-1) = 2*K(m)*K(m-1)/(K(m)+K(m-1))/Dx;
b(m) = -N(m)*Dx-2*phi0R*K(m)/Dx;
%Mass balance for i-th intermediate block:
for k=2:m-1
    A(k,k) = -2*(K(k)*K(k-1)/(K(k)+K(k-1))+K(k)*K(k+1)/(K(k)+K(k+1)))/Dx; 
    A(k,k-1) = 2*K(k)*K(k-1)/(K(k)+K(k-1))/Dx;
    A(k,k+1) = 2*K(k)*K(k+1)/(K(k)+K(k+1))/Dx;
    b(k) = -N(k)*Dx;
end