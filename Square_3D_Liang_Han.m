

%% RBF-FD : NSCH-Coupled
clear all
clc
close all
format shorte
%% initial points
Ni = 20;
x_min = -0.3;
x_max = 0.3;
y_min = -0.3;
y_max = 0.3;
z_min = -0.3;
z_max = 0.3;

h1 = (x_max-x_min)/Ni;
h2 = (y_max-y_min)/Ni;
h3 = (z_max-z_min)/Ni;

x = x_min:h1:x_max;
y = y_min:h2:y_max;
z = z_min:h3:z_max;
[X,Y,Z] = meshgrid(x,y,z);
x = X(:);   y = Y(:);  z = Z(:);
M = length(x);
%%
N = length(x);
Boundary=find(x==x_min | x==x_max | y==x_min | y==x_max| z==z_min | z==z_max);
%%
west=find(y==y_min & x~=x_min & x~=x_max & z~=z_min & z~=z_max);
east=find(y==y_max & x~=x_min & x~=x_max & z~=z_min & z~=z_max);
south=find(z==z_min & x~=x_min & x~=x_max & y~=y_min & y~=y_max);
north=find(z==z_max & x~=x_min & x~=x_max & y~=y_min & y~=y_max);
front=find(x==x_max & y~=y_min & y~=y_max & z~=z_min & z~=z_max);
back=find(x==x_min & y~=y_min & y~=y_max & z~=z_min & z~=z_max);

%%
% border1=[west;east;south;north];
% corner=(1:N)';
% corner([border1;Interior])=[];
%%
eta = 15e-2;       % viscosity parameter
rho = 1;       % volume-averaged density
gamma = 1e-2;     % surface tension between the two fluid phases
eps = 1e-2;    % dimensionless parameter
Mob = 1e-2;       % Mobility Function
%% Time discrete scheme
dt = 1e-1;
T_final = 1;
T=0:dt:T_final;
Nt = length(T);
TSCREEN = 1;
%% Basis Functions and Required Matrices
phi = @(r,c) exp(-(c*r).^2);
drphi = @(r,rx,c) -2*c^2*rx.*exp(-(c*r).^2);
d2rphi = @(r,rx,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rx.^2));

%%

ns = 10;
% ci = zeros(M,1);
F = zeros(M,ns);
C = zeros(M,ns);
for i=1:M
    x_center = x(i);
    y_center = y(i);
    z_center = z(i);
    
    rd = sqrt((x_center-x).^2+(y_center-y).^2+(z_center-z).^2);
    [rd,ix] = sort(rd);
    F(i,:) = rd(1:ns);
    C(i,:) = ix(1:ns);
%     ci(i) = rd(ns)*sqrt(ns)*rd(2)/(0.02*(sum(F(i,1:ns))));
end

%%
c = 10;
minK = 1e10;
maxK = 1e12;
dc = 0.1;

%%
wxx = zeros(M);
wyy = zeros(M);
wzz = zeros(M);
wx = zeros(M);
wy = zeros(M);
wz = zeros(M);


for i=1:M
    pn = C(i,:);
    rx = x(pn)-x(pn)';
    ry = y(pn)-y(pn)';
    rz = z(pn)-z(pn)';
    
    r = sqrt(rx.^2+ry.^2+rz.^2);
    
%      A_local = phi(r,c(i));
     D_local = sqrt(rx.^2+ry.^2+rz.^2);
     K = 1;
    while (K<minK || K>maxK)                           
        A_local=phi(D_local,c);
        [~,Ze,~]=svd(A_local);
        K = Ze(1,1)/Ze(ns,ns);
        if K<minK  
            c = c - dc;
        elseif K>maxK 
            c = c + dc;
        end
    end

    
        
    Bx = (drphi(sqrt((x(i)-x(pn)).^2+(y(i)-y(pn)).^2)+(z(i)-z(pn)).^2,x(i)-x(pn),c));
    By = (drphi(sqrt((x(i)-x(pn)).^2+(y(i)-y(pn)).^2)+(z(i)-z(pn)).^2,y(i)-y(pn),c));
    Bz = (drphi(sqrt((x(i)-x(pn)).^2+(y(i)-y(pn)).^2)+(z(i)-z(pn)).^2,z(i)-z(pn),c));
    Bxx = (d2rphi(sqrt((x(i)-x(pn)).^2+(y(i)-y(pn)).^2)+(z(i)-z(pn)).^2,x(i)-x(pn),c));
    Byy = (d2rphi(sqrt((x(i)-x(pn)).^2+(y(i)-y(pn)).^2)+(z(i)-z(pn)).^2,y(i)-y(pn),c));
    Bzz = (d2rphi(sqrt((x(i)-x(pn)).^2+(y(i)-y(pn)).^2)+(z(i)-z(pn)).^2,z(i)-z(pn),c));
    
    
    wx(i,pn) = pinv(A_local)*Bx;
    wy(i,pn) = pinv(A_local)*By;
    wz(i,pn) = pinv(A_local)*Bz;
    wxx(i,pn) = pinv(A_local)*Bxx;
    wyy(i,pn) = pinv(A_local)*Byy;
    wzz(i,pn) = pinv(A_local)*Bzz;

 
end

%%
Dx = wx;  Dy = wy;  Dz = wz;
Dxx = wxx;  Dyy = wyy;  Dzz = wzz;
%%
U = zeros(N,1);
V = zeros(N,1);
W = zeros(N,1);

I = eye(N);
O = zeros(N);
%% Pressure Equation
DD=Dxx+Dyy+Dzz;

Lpsi=Dxx+Dyy+Dzz;
Lpsi(north,:)=Dy(north,:);
Lpsi(south,:)=Dy(south,:);
Lpsi(east,:)=Dx(east,:);
Lpsi(west,:)=Dx(west,:);
Lpsi(front,:)=Dz(front,:);
Lpsi(back,:)=Dx(back,:);

%% Initial zero point Functions


Cn = zeros(N,1);
% Cn = tanh((0.2-sqrt((x-0.5).^2+(y-2.4).^2))/(sqrt(2)*eps))+tanh((2-y)/(sqrt(2)*eps));
% Cn = tanh(2*(0.25-sqrt((x-0.5).^2+(y-0.5).^2))/(sqrt(2)*eps));
% Cn = tanh(2*(0.25-sqrt((x-0.5-0.25).^2+(y-0.5).^2))/eps) + tanh(2*(0.25-sqrt((x-0.5+0.25).^2+(y-0.5).^2))/eps);


% Cn = tanh(2*(0.3-abs(x+y+z-1)-abs(x-y-z))/eps);




Cn = -tanh((1/2 *(abs(x)+abs(y)+abs(z))-0.1)/(sqrt(2)*eps));
% Cn = tanh(2*(0.3-abs(x+y+z-1)-abs(x-y)-abs(x-z))/eps);


isosurface(X,Y,Z,reshape(Cn,size(X)))


%%
Fef0 = (1/4)*(1-Cn.^2).^2; %quartic homogeneous free energy density function
    Fef0_Cn = Cn.^3 - Cn;     % Derivate to Cn
    mu = -gamma*eps*DD*Cn + (gamma/3)*(Cn.^3 - Cn);
%%    
    A1 = I ;
    A2 = -dt*DD;
    A3 = gamma*eps*DD;
    A4 = I;
 
    LCH = [A1,A2 ; A3,A4] ; %left side of Cahn-Hilliard system
    
    LCH([east;west],1:N)=Dx([east;west],:);
    LCH([south;north],1:N)=Dy([south;north],:); 
    LCH([front;back],1:N)=Dz([front;back],:);
    
    LCH([south;north],N+1:end)=O([south;north],:);
    LCH([east;west],N+1:end)=O([east;west],:);
    LCH([front;back],N+1:end)=O([front;back],:);
    
    LCH([east;west]+N,N+1:end)=Dx([east;west],:);
    LCH([south;north]+N,N+1:end)=Dy([south;north],:);
    LCH([front;back]+N,N+1:end)=Dz([front;back],:);
    
    LCH([east;west]+N,1:N)=O([east;west],:);
    LCH([south;north]+N,1:N)=O([south;north],:);
    LCH([front;back]+N,1:N)=O([front;back],:);
    
    %%
k=0;

for i=1:Nt-1
    k=k+1;
    
    
    %%
    U_star = U - dt*U.*(Dx*U) - dt*V.*(Dy*U) - dt*W.*(Dz*U) + eta*(dt/rho)*DD*U - (dt/rho)*Cn.*(Dx*mu);
    
    V_star = V - dt*U.*(Dx*V) - dt*V.*(Dy*V) - dt*W.*(Dz*V) + eta*(dt/rho)*DD*V - (dt/rho)*Cn.*(Dy*mu);
    
    W_star = W - dt*U.*(Dx*W) - dt*V.*(Dy*W) - dt*W.*(Dz*W) + eta*(dt/rho)*DD*W - (dt/rho)*Cn.*(Dz*mu);
    
    
    U_star(Boundary)=0;
    
    V_star(Boundary)=0;
    
    W_star(Boundary)=0;
    
    %%
    RHS_Poisson = (rho/dt)*(Dx*U_star + Dy*V_star + Dz*W_star);
    
    RHS_Poisson(Boundary) = 0;
    
    P = Lpsi\RHS_Poisson;
    %%
    U = U_star - (dt/rho)*Dx*P;
    
    V = V_star - (dt/rho)*Dy*P;
    
    W = W_star - (dt/rho)*Dz*P;
    
    U(Boundary)=0;
    
    V(Boundary)=0;
    
    W(Boundary)=0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    RCH_up = Cn - dt*Dx*(Cn.*U) - dt*Dy*(Cn.*V) - dt*Dz*(Cn.*W); 
    RCH_down = (gamma/3).*(Cn.^3 - Cn);
    
    RCH = [RCH_up;RCH_down];
    
    CNM = LCH\RCH;
    
    
    Cn = CNM(1:N);  %update Cn
    mu = CNM(N+1:end);
    
     
    if (k==TSCREEN)
            

        isosurface(X,Y,Z,reshape(Cn,size(X)))

        colormap(gca,"jet")
        colorbar


        drawnow
        k=0;
%         pause
    end
    
    T(i+1)
end
