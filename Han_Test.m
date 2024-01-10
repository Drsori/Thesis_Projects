

%% RBF-FD : NSCH-Coupled
clear all
clc
close all
format shorte
%% initial points
Ni = 32;
x_min = 0;
x_max = 1;
y_min = 0;
y_max = 1;
h = (x_max-x_min)/Ni;
z = (y_max-y_min)/Ni;
x = x_min:h:x_max;
y = y_min:z:y_max;
[X,Y] = meshgrid(x,y);
x = X(:);   y = Y(:);
M = length(x);
%%
N = length(x);
Boundary=find(x==x_min | x==x_max | y==x_min | y==x_max);
Interior=find(x~=x_min | x~=x_max | y~=x_min | y~=x_max);
%%
Boundary1=find(x==x_min);
Boundary2=find(x==x_max);
Boundary3=find(y==y_min);
Boundary4=find(y==y_max);
%%
west=find(x==x_min & y~=x_min & y~=x_max);
east=find(x==x_max & y~=x_min & y~=x_max);
south=find(y==x_min & x~=x_min & x~=x_max);
north=find(y==x_max & x~=x_min & x~=x_max);
%%
% border1=[west;east;south;north];
% corner=(1:N)';
% corner([border1;Interior])=[];
%%
Re = 10;
We = 200;     % Weber number
eps = 5e-3;    % dimensionless parameter
%% Time discrete scheme
dt = 5e-3;
T_final = 10;
T=0:dt:T_final;
Nt = length(T);
TSCREEN = 1;
%% Basis Functions and Required Matrices
phi = @(r,c) exp(-(c*r).^2);
drphi = @(r,rx,c) -2*c^2*rx.*exp(-(c*r).^2);
d2rphi = @(r,rx,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rx.^2));

% phi = @(r,c) sqrt(r.^2+c.^2);
% drphi = @(r,rh,s) rh./sqrt(r.^2+s.^2);
% d2rphi = @(r,rh,s) 1./sqrt(r.^2+s.^2)-(rh.^2)./(r.^2+s.^2).^(1.5);


%%

ns = 8;
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
    rd = sqrt(rx.^2+ry.^2);
    
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
A = phi(rd,c);
Dx = wx;  Dy = wy;
Dxx = wxx;  Dyy = wyy;
%%
DD = Dxx+Dyy;

%% Pressure Equation
Lpsi=Dxx+Dyy;
Lpsi(north,:)=Dy(north,:);
Lpsi(south,:)=Dy(south,:);
Lpsi(east,:)=Dx(east,:);
Lpsi(west,:)=Dx(west,:);
%%
k=0;

%% Initial zero point Functions
Cn = zeros(N,1);
% Cn = -tanh((abs(x+y-1)+abs(x-y)-0.2)/(sqrt(2)*eps));
Cn = (625e-4)*(1+tanh(3*x/eps)).*(1+tanh(3*(x-1)/eps)).*(1+tanh(3*y/eps)).*(1+tanh(3*(y-1)/eps));
%%%%%%%
%% random number between a,b :  r = a + (b-a)*rand(N,1);

% r = -5e-2 + (5e-2+5e-2)*rand(N,1);
% Cn0 = @(x,y,T) -5e-2 + r ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = zeros(N,1);
V = zeros(N,1);
I = eye(N);
Ze =zeros(N);

% U(Boundary4)=x(Boundary4).*(1-x(Boundary4));
% V(Boundary4)=y(Boundary4).*(1-y(Boundary4));

% mesh(X,Y,reshape(Cn,size(X))); shading interp, drawnow
 
%%
    Fef0 = (1/4)*(1-Cn.^2).^2; %quartic homogeneous free energy density function
    Fef0_Cn = Cn.^3 - Cn;     % Derivate to Cn
    mu = Fef0_Cn - (eps^2)*DD*Cn;
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Mobility Update
      Mob = (0.1)*sqrt((1-Cn.^2).^2 +(eps^2));     % Mobility Function
%% Chemical Potential and Cahn Hilliard System
    A1 = I ;
    A2 = -dt*Mob.*DD;
    A3 = ((eps)^2)*DD;
    A4 = I;
 
    LCH = [A1,A2 ; A3,A4] ; %left side of Cahn-Hilliard system
    
    LCH([south;north],1:N)=Dy([south;north],:);
    LCH([east;west],1:N)=Dx([east;west],:);
    LCH([south;north],1+N:end)=Ze([south;north],:);
    LCH([east;west],1+N:end)=Ze([east;west],:);

    LCH([south;north]+N,1+N:end)=Dy([south;north],:);
    LCH([south;north]+N,1:N)=Ze([south;north],:);
    LCH([east;west]+N,1+N:end)=Dx([east;west],:);
    LCH([east;west]+N,1:N)=Ze([east;west],:);


    

for i=1:Nt-1
    i*dt
    k=k+1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U_star = U - dt*U.*(Dx*U) - dt*V.*(Dy*U) + (dt/Re)*(Dxx + Dyy)*U - (dt/(eps*We))*Cn.*(Dx*mu);
    
    V_star = V - dt*U.*(Dx*V) - dt*V.*(Dy*V) + (dt/Re)*(Dxx + Dyy)*V - (dt/(eps*We))*Cn.*(Dy*mu);
    
    
    U_star(Boundary)=0;
    
    V_star(Boundary)=0;
    
    
    %%
    RHS_Poisson = (1/dt)*(Dx*U_star + Dy*V_star);
    
    RHS_Poisson(Boundary) = 0;
    
    P = Lpsi\RHS_Poisson;
    %%
    U = U_star - dt*Dx*P;
    
    V = V_star - dt*Dy*P;
    
    U(Boundary)=0;
    
    V(Boundary)=0;
    
    %%
    RCH_up = Cn - dt*Dx*(Cn.*U) + dt*Dy*(Cn.*V);
    RCH_Down = Cn.^3 - Cn;

    RCH = [RCH_up;RCH_Down];
    
    CNM = LCH\RCH;
    
    
    Cn = CNM(1:N);  %update Cn
    mu = CNM(N+1:end);
 
    
%     mesh(X,Y,reshape(Cn,size(X))); shading interp, drawnow
    
       if (k==TSCREEN)
        %        Go back in real space omega in real space for plotting
%         pcolor(X,Y,reshape(Phi,size(X))); shading interp

%             colormap(jet)
%             h=figure;
            contourf(X,Y,reshape(Cn,size(X)));
            colormap(jet)
            colorbar
            drawnow
%             saveas(h,[pwd sprintf( './figs_hen/FIG%d_(t=%f).png',i,T(i+1))]);
%             close(h)


        %         plot(0)
        %         hh=streamslice(X,Y,reshape(U,size(X)),reshape(V,size(X)),25,'k','noarrow');
        %         set( hh, 'Color', 'b');
        %         hold off
        %     drawnow
            k=0;
%         pause
        end
    
    T(i+1)
end
