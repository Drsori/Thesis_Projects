

%% RBF-FD : NSCH-Coupled
clear all
clc
close all
format shorte
%% initial points
Ni = 40;
x_min = 0;
x_max = 2;
y_min = 0;
y_max = 2;
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
eta = 1e-1;       % viscosity parameter
rho = 1;       % volume-averaged density
gamma = 1e-1;     % surface tension between the two fluid phases
eps = 1e-2;    % dimensionless parameter
Mob = 1e-1;       % Mobility Function
%% Time discrete scheme
dt = 1e-2;
T_final = 10;
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


% Cn = zeros(N,1);
% teta = zeros(N,1);
% for i=1:N
%    teta(i) = atan(y(i)/x(i));
% end
% teta(1)=teta(2);
% R = sqrt(x.^2 + y.^2);
% Cn = tanh(1/(sqrt(2)*eps) *(R - 0.25 - (1+cos(4*teta))/16 ));


% Cn = 1+sin((0.28-sqrt((x.^2-0.7).^2+(y.^2-1).^2))/(sqrt(2)*eps));


%% RING
Cn = tanh(1/(sqrt(2)*eps) *((x.^2)/(eps) + (y.^2)/(eps) -1).*((x.^2)/(eps) + (y.^2)/(eps) -1));

% mesh(X,Y,reshape(Cn,size(X))); shading interp, drawnow



% rnd = -1 + (2)*rand(N,1);
% Cn = 0.9*(y/Ly - 0.5) + (1e-3)*rnd;
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
    LCH([south;north],N+1:end)=Z([south;north],:);
    LCH([east;west],N+1:end)=Z([east;west],:);
    
    LCH([east;west]+N,N+1:end)=Dx([east;west],:);
    LCH([south;north]+N,N+1:end)=Dy([south;north],:);
    LCH([east;west]+N,1:N)=Z([east;west],:);
    LCH([south;north]+N,1:N)=Z([south;north],:);
%%
k=0;

for i=1:Nt-1
    k=k+1;
    
    
    %%
    U_star = U - dt*U.*(Dx*U) - dt*V.*(Dy*U) + eta*(dt/rho)*DD*U - (dt/rho)*Cn.*(Dx*mu);
    
    V_star = V - dt*U.*(Dx*V) - dt*V.*(Dy*V) + eta*(dt/rho)*DD*V - (dt/rho)*Cn.*(Dy*mu);
    
    U_star(Boundary)=0;
    
    V_star(Boundary)=0;
    
    U_star(north)=1;
    
    %%
    RHS_Poisson = (rho/dt)*(Dx*U_star + Dy*V_star);
    
    RHS_Poisson(Boundary) = 0;
    
    P = Lpsi\RHS_Poisson;
    %%
    U = U_star - (dt/rho)*Dx*P;
    
    V = V_star - (dt/rho)*Dy*P;
    
    U(Boundary)=0;
    
    V(Boundary)=0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    RCH_up = Cn - dt*Dx*(Cn.*U) + dt*Dy*(Cn.*V); 
    RCH_down = (gamma/3).*(Cn.^3 - Cn);
    
    RCH = [RCH_up;RCH_down];
    
    CNM = LCH\RCH;
    
    
    Cn = CNM(1:N);  %update Cn
    mu = CNM(N+1:end);
    
     
    if (k==TSCREEN)
        %        Go back in real space omega in real space for plotting
%         pcolor(X,Y,reshape(Cn,size(X))); shading interp
            colormap(jet)
%             h=figure;
            contourf(X,Y,reshape(Cn,size(X)),1);
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
