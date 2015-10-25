function K0=q2()
    close all;
    m = 200;
    K0 = zeros(m,1);
    K0(1:20,1)=1e-6;
    K0(21:40,1)=1e-6;
    K0(41:60,1)=1e-6;
    K0(61:80,1)=2e-6;
    K0(81:100,1)=2e-6;
    K0(101:120,1)=2e-6;
    K0(121:140,1)=3e-6;
    K0(141:160,1)=2e-6;
    K0(161:180,1)=2e-6;
    K0(181:200,1)=3e-6;
    Dx = 1/m;
    phi0L=1; phi0R=1;
    N = 1E-5*ones(m,1);
    [A,b] = makeab1(K0,N,Dx,phi0L,phi0R);
    phi = A\b;
    x = linspace(0,1,m);
    data = load('head_obs.txt');
    x_obs = data(:,1); phi_obs = data(:,2);
    figure;
    plot(x,phi,'k.-');
    hold on;
    plot(x_obs,phi_obs,'ro');
    

end