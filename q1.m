function q1()
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





end