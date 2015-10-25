function q1()
    clc; clear;
    close all;
    m = 200;
    K = 1e-6*ones(m,1);
    Dx = 1/m;
    phi0L=1; phi0R=1;
    N = 1E-5*ones(m,1);
    [A,b] = makeab1(K,N,Dx,phi0L,phi0R);
    phi = A\b;
    x = linspace(0,1,m);
    data = load('head_obs.txt');
    x_obs = data(:,1); phi_obs = data(:,2);
    figure;
    plot(x,phi,'k.-');
    hold on;
    plot(x_obs,phi_obs,'ro');

    %% Verify if cderab1.m work as intended
    ind = 1;
    [Ak_analytical,bk_analytical] = cderab1(ind,K,N,Dx,phi0L,phi0R);
    dk = 1e-12; dK = zeros(m,1); dK(ind) = dk;
    [A_dk,b_dk] = makeab1(K+dK,N,Dx,phi0L,phi0R);
    Ak_numerical = (A_dk - A )/dk;
    bk_numerical = (b_dk - b )/dk;
    normest(Ak_numerical-Ak_analytical)/normest(Ak_numerical)
    normest(bk_numerical-bk_analytical)/normest(bk_numerical)
end