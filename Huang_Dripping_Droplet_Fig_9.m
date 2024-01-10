

%% RBF-FD : NSCH-Coupled
clear all
clc
close all
format shorte
%% initial points
Ni = 20;
x_min = 0;
x_max = 1;
y_min = 0;
y_max = 4;
h = (x_max-x_min)/Ni;
z = (x_max-x_min)/Ni;
x = x_min:h:x_max;
y = y_min:z:y_max;
[X,Y] = meshgrid(x,y);
x = X(:);   y = Y(:);
M = length(x);
%%
N = length(x);
Boundary=find(x==x_min | x==x_max | y==y_min | y==y_max);
Interior=find(x~=x_min | x~=x_max | y~=x_min | y~=x_max);
%%
Boundary1=find(x==x_min);
Boundary2=find(x==x_max);
Boundary3=find(y==y_min);
Boundary4=find(y==y_max);
%%
west=find(x==x_min & y~=y_min & y~=y_max);
east=find(x==x_max & y~=y_min & y~=y_max);
south=find(y==y_min & x~=x_min & x~=x_max);
north=find(y==y_max & x~=x_min & x~=x_max);
%%
% border1=[west;east;south;north];
% corner=(1:N)';
% corner([border1;Interior])=[];
%%
rho1 = 1.6;
rho2 = 1;
rho_star = rho2;      
Re = 30;       
We = 100;       
eps = 15e-3;
Pe = 1/eps;
%% Time discrete scheme
dt = 1e-2;
T_final = 20;
T=0:dt:T_final;
Nt = length(T);
TSCREEN = 1;
%% Basis Functions and Required Matrices
phi = @(r,c) exp(-(c*r).^2);
drphi = @(r,rx,c) -2*c^2*rx.*exp(-(c*r).^2);
d2rphi = @(r,rx,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rx.^2));

%%

ns = 15;
% ci = zeros(M,1);
F = zeros(M,ns);
C = zeros(M,ns);
for i=1:M
    x_center = x(i);
    y_center = y(i);
    rd = sqrt((x_center-x).^2+(y_center-y).^2);
    [rd,ix] = sort(rd);
    F(i,:) = rd(1:ns);
    C(i,:) = ix(1:ns);
%     ci(i) = rd(ns)*sqrt(ns)*rd(2)/(0.02*(sum(F(i,1:ns))));
end

%%
c = 10;
minK = 1e10;
maxK = 1e14;
dc = 0.1;

%%
wxx = zeros(M);
wyy = zeros(M);
wx = zeros(M);
wy = zeros(M);


for i=1:M
    pn = C(i,:);
    rx = x(pn)-x(pn)';
    ry = y(pn)-y(pn)';
    r = sqrt(rx.^2+ry.^2);
    
%      A_local = phi(r,c(i));
     D_local = sqrt(rx.^2+ry.^2);
     K = 1;
    while (K<minK || K>maxK)                           
        A_local=phi(D_local,c);
        [~,Z,~]=svd(A_local);
        K = Z(1,1)/Z(ns,ns);
        if K<minK  
            c = c - dc;
        elseif K>maxK 
            c = c + dc;
        end
    end

    
        
    Bx = (drphi(sqrt((x(i)-x(pn)).^2+(y(i)-y(pn)).^2),x(i)-x(pn),c));
    By = (drphi(sqrt((x(i)-x(pn)).^2+(y(i)-y(pn)).^2),y(i)-y(pn),c));
    Bxx = (d2rphi(sqrt((x(i)-x(pn)).^2+(y(i)-y(pn)).^2),x(i)-x(pn),c));
    Byy = (d2rphi(sqrt((x(i)-x(pn)).^2+(y(i)-y(pn)).^2),y(i)-y(pn),c));
    
    wx(i,pn) = pinv(A_local)*Bx;
    wy(i,pn) = pinv(A_local)*By;
    wxx(i,pn) = pinv(A_local)*Bxx;
    wyy(i,pn) = pinv(A_local)*Byy;
    

 
end

%%
Dx = wx;  Dy = wy;
Dxx = wxx;  Dyy = wyy;
%%
U = zeros(N,1);
V = zeros(N,1);
I = eye(N);
Z = zeros(N);
%% Pressure Equation
DD=Dxx+Dyy;

Lpsi=Dxx+Dyy;
Lpsi(north,:)=Dy(north,:);
Lpsi(south,:)=Dy(south,:);
Lpsi(east,:)=Dx(east,:);
Lpsi(west,:)=Dx(west,:);

%% Initial zero point Functions


Cn = zeros(N,1);

R = sqrt((x-0.5).^2 + (y-4.2).^2);
Cn = tanh((0.5- R)/(eps*sqrt(2)));

% pcolor(X,Y,reshape(Cn,size(X))); shading interp
% axis equal
% colormap(jet)
%%
%%
    Fef0 = (1/4)*(1-Cn.^2).^2; %quartic homogeneous free energy density function
    Fef0_Cn = Cn.^3 - Cn;     % Derivate to Cn
    Mob = (1-Cn.^2).^2 ;     % Mobility Function
    mu = Fef0_Cn - (eps^2)*(Dxx+Dyy)*Cn;
%%    
    A1 = I ;
    A2 = -(dt/Pe)*DD;
    A3 = (eps^2)*DD;
    A4 = I;
    
    A1(east,:)  = -Dx(east,:);
    A1(west,:)  =  Dx(west,:);
    A1(south,:) = -Dy(south,:);
    A1(north,:) =  Dy(north,:);

    A4(east,:)  = -Dx(east,:);
    A4(west,:)  =  Dx(west,:);
    A4(south,:) = -Dy(south,:);
    A4(north,:) =  Dy(north,:);

    LCH = [A1,A2 ; A3,A4] ; %left side of Cahn-Hilliard system
    
    
    
k=0;
for i=1:Nt-1
    k=k+1;
    
    rho = (rho1/2)*(1 + Cn) + (rho2/2)*(1 + Cn);
    
    %%
    U_star = U - dt*U.*(Dx*U) - dt*V.*(Dy*U) + (dt/(rho_star*Re))*DD*U - (dt/(rho_star*We))*Cn.*(Dx*mu);
    
    V_star = V - dt*U.*(Dx*V) - dt*V.*(Dy*V) + (dt/(rho_star*Re))*DD*V - (dt/(rho_star*We))*Cn.*(Dy*mu) - dt*(rho/rho_star - 1);
    
    U_star(Boundary)=0;
    
    V_star(Boundary)=0;
    
    
    %%
    RHS_Poisson = (rho_star/dt)*(Dx*U_star + Dy*V_star);
    
    RHS_Poisson(Boundary) = 0;
    
    P = Lpsi\RHS_Poisson;
    %%
    U = U_star - (dt/rho)*Dx*P;
    
    V = V_star - (dt/rho)*Dy*P;
    
    U(Boundary)=0;
    
    V(Boundary)=0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    RCH_up = Cn - dt*Dx*(Cn.*U) + dt*Dy*(Cn.*V); 
    RCH_down = (Cn.^3 - Cn);
    
    
    RCH = [RCH_up;RCH_down];
    
    CNM = LCH\RCH;
    
    
    Cn = CNM(1:N);  %update Cn
    mu = CNM(N+1:end);
    
     
    if (k==TSCREEN)
        %        Go back in real space omega in real space for plotting
            colormap(jet)        
            pcolor(X,Y,reshape(Cn,size(X))); shading interp
            axis equal
%             h=figure;
%             contourf(X,Y,reshape(Cn,size(X)),4);
            drawnow
%             saveas(h,[pwd sprintf( './figs/FIG%d_(t=%f).png',i,T(i+1))]);
%             close(h)
        %         plot(0)
        %         hh=streamslice(X,Y,reshape(U,size(X)),reshape(V,size(X)),25,'k','noarrow');
        %         set( hh, 'Color', 'b');
        %         hold off
        drawnow
        k=0;
%         pause
    end
    
    T(i+1)
end

