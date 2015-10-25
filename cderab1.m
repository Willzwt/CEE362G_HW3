function [Ak,bk] = cderab1(k,K,N,Dx,phi0L,phi0R)
%[Ak,bk] = cderab1(k,K,N,Dx,phi0L,phi0R)
%computes derivative of A and b with respect to K(k)
%  (see also makeab1.m)
%where:
%  k is integer 
%  K is m by 1 vector of conductances
%  N is m by 1 vector of recharge rates (included for generality)
%  Dx is size of blocks (same for all)
%  phi0L and phi0R are left and right boundary head values.  
%    
%Note that Ak and bk are sparse.  

%pkk  May 5, 2009

m = length(K);
Ak = sparse(zeros(m,m)); bk = sparse(m,1);
if k==1
    % For first-block conductivity:
    Ak(1,1) = -2*(1+K(2)^2/(K(1)+K(2))^2)/Dx;
    Ak(1,2) = 2*K(2)^2/(K(1)+K(2))^2/Dx;
    Ak(2,1) = Ak(1,2);
    Ak(2,2) = -2*K(2)^2/(K(2)+K(1))^2/Dx; 
    bk(1) = -2*phi0L/Dx;
elseif k==m
    %Last-block conductivity:
    Ak(m,m) = -2*(1+K(m-1)^2/(K(m)+K(m-1))^2)/Dx;
    Ak(m,m-1) = 2*K(m-1)^2/(K(m)+K(m-1))^2/Dx;
    Ak(m-1,m) = Ak(m,m-1);
    Ak(m-1,m-1) = -2*K(m-1)^2/(K(m)+K(m-1))^2/Dx; 
    bk(m) = -2*phi0R/Dx;
elseif k>1 & k<m
    %For k-th intermediate block:
    Ak(k,k) = -2*(K(k-1)^2/(K(k)+K(k-1))^2+K(k+1)^2/(K(k)+K(k+1))^2)/Dx;
    Ak(k-1,k-1) = -2*K(k-1)^2/(K(k)+K(k-1))^2/Dx;
    Ak(k+1,k+1) = -2*K(k+1)^2/(K(k)+K(k+1))^2/Dx;
    Ak(k,k-1) = 2*K(k-1)^2/(K(k)+K(k-1))^2/Dx;
    Ak(k-1,k) = Ak(k,k-1);
    Ak(k,k+1) = 2*K(k+1)^2/(K(k)+K(k+1))^2/Dx;
    Ak(k+1,k) = Ak(k,k+1);
end