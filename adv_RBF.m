

%% Advection RBF-FD
clear all
clc
format shorte
N = 20;
a = 0;   b = 1;
eps = 5e-5;
h = (b-a)/N;
T = 1;
Nt = 1e5;
dt = T/Nt;
x = a:h:b;
x = x';
M = length(x);

%%
ns = 5;
c = zeros(M,1);
F = zeros(M,ns);
C = zeros(M,ns);
%% Basis Functions and Required Matrices
phi = @(r,c) exp(-c^2.*r.^2);
dxphi = @(r,rx,c) -2*c^2*rx.*exp(-c^2.*r.^2);
dxxphi = @(r,rx,c) (-2*c^2+4*c^4.*rx.^2).*exp(-c^2.*r.^2);
%%
for i=1:M
    x_center = x(i);
    rd = sqrt((x_center-x(:)).^2);
    [rd,ix] = sort(rd);
    F(i,:) = rd(1:ns);
    C(i,:) = ix(1:ns);
    c(i) = rd(ns)*sqrt(ns)*rd(2)/(0.02*(sum(F(i,1:ns))));
end

%%
I = ones(1,ns);
wx1 = zeros(M);
wx2 = zeros(M);

%%
for i=1:M
    pn = C(i,:);
    rx = x(pn)*I-(x(pn)*I)';
    r = sqrt(rx.^2);
    A_local = phi(r,c(i));
    B_local_1 = dxphi(sqrt((x(i)-x(pn)).^2),x(i)-x(pn),c(i));
    B_local_2 = dxxphi(sqrt((x(i)-x(pn)).^2),x(i)-x(pn),c(i));
    wx1(i,pn) = pinv(A_local)*B_local_1;
    wx2(i,pn) = pinv(A_local)*B_local_2;
end

%%
u = sin(2*pi.*x)+0.5.*sin(pi.*x);
RHS = @(u) -u.*(wx1*u)+eps*wx2*u;
%% RK4
for k=1:Nt
    t = k*dt
    k1 = dt*RHS(u); 
    k2 = dt*RHS(u+k1/2);
    k3 = dt*RHS(u+k2/2);
    k4 = dt*RHS(u+k3);
    u = u + (1/6)*(k1+2*(k2+k3)+k4);
    u(1) = 0;
    u(end) = 0;
    if mod(k,500)==0
        plot(x,u','.-'),drawnow
    end
end
plot(x,u')


