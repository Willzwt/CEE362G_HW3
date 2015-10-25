function q3()
close all; clc; clear;
m = 200;
Dx = 1/m;
phi0L=1; phi0R=1;
N = 1E-5*ones(m,1);
x = linspace(0,1,m)';
%% Read observation data
data = load('head_obs.txt');
x_obs = data(:,1); phi_obs = data(:,2);
n = length(phi_obs);
ind = dsearchn(x,x_obs);
%% Function to calculate the Jocabian matrix H given K
    function [H,phi] = getJacobian(K)
        [A,b] = makeab1(K,N,Dx,phi0L,phi0R);
        phi = A\b;
        H = zeros(n,m);
        I = eye(m);
        lamba = (A')\I(:,ind);
        for row = 1:n
            % Fill the Jacobian matrix row by row
            for col = 1:m
                [Ak,bk] = cderab1(col,K,N,Dx,phi0L,phi0R);
                As = K(col)*Ak;
                bs = K(col)*bk;
                H(row,col) = lamba(:,row)' * (bs-As*phi);
            end
        end
    end
%% Cross validate function
    function [] = crossValidate(K)
        [A,b] = makeab1(K,N,Dx,phi0L,phi0R);
        phi_temp = A\b;
        plot(x,phi_temp,'k.-');
        hold on;
        plot(x_obs,phi_obs,'ro');
    end
%% Prior knowledge of model K
K0 = 2e-6*ones(m,1);
s0 = log(K0);
std_s = 0.1*mean(s0);var_s = std_s^2;
l_s = 1/5;
R = 0.0005*mean(phi_obs)*eye(n);
X = ones(m,1);
%% Solving the problem iteratively
Q = gaussianQ(m,x,var_s,l_s);
[H,phi] = getJacobian(K0);
y = phi_obs - phi(ind) + H*s0;
[s,V,LAMDA,MU] = GenLinInv(y,H,R,X,Q);
K = exp(s);
crossValidate(K);
eps = 1e-3;
% iterate until convergence
figure(1);
crossValidate(K0);
hold on;
figure(2);
hold on;
plot(x,s0);
n_iter = 0;
while (norm(K-K0)/norm(K0)>eps)
    K0 = K; s0 = s;
    [H,phi] = getJacobian(K);
    y = phi_obs - phi(ind) + H*s0;
    [s,V,LAMDA,MU] = GenLinInv(y,H,R,X,Q);
    K = exp(s);
    figure(1);
    crossValidate(K);
    figure(2);
    plot(x,s);
    n_iter = n_iter + 1
end
    

    


end



function Q = ExponentialCov(n,t,v,l)

h = abs(repmat(t,1,n) - repmat(t',n,1));
Q = v*exp(-h/l);

end


function Q=linearQ(n,x,var,l)

theta = var/l;
h=abs(repmat(x,1,n)-repmat(x',n,1));
Q=-theta * h;

end

function Q=gaussianQ(n,x,var,L)

h=abs(repmat(x,1,n)-repmat(x',n,1));
Q=var*exp(-(h.^2/L^2));

end

function Q=smoothCov(n,t)

h = abs(repmat(t,1,n) - repmat(t',n,1));
Q = h.^3;

end