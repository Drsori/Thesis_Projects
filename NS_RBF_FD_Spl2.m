clear all
clc
close all
format shorte
%% Collocation Points
x_min = 0;
x_max = 1;
y_min = 0;
y_max = 1;
Ni = 20;
h = (x_max-x_min)/Ni;
x = x_min:h:x_max;
y = y_min:h:y_max;
[X,Y] = meshgrid(x,y);
x = X(:);  y = Y(:);

Boundary = find( x==x_min | x==x_max | y==y_min | y==y_max );
Interior = find( x~=x_min | x~=x_max | y~=y_min | y~=y_max );
N = length(x);
[M,~] = size(x);
%%
west=find(x==x_min & y~=y_min & y~=y_max);
east=find(x==x_max & y~=y_min & y~=y_max);
south=find(y==y_min & x~=x_min & x~=x_max);
north=find(y==y_max & x~=x_min & x~=x_max);
%%
border1=[west;east;south;north];
corner=(1:N)';
corner([border1;Interior])=[];
%% Time discrete scheme
dt=1e-4;
T_final=10;
T=0:dt:T_final;
Nt=length(T);
TSCREEN=10;
%% Basis Functions and Required Matrices
Phi    = @(r,c) exp(-(c*r).^2);
drphi  = @(r,rh,c) -2*c^2*rh.*exp(-(c*r).^2);
d2rphi = @(r,rh,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rh.^2));


%%
ns = 5;
ci = zeros(M,1);
F = zeros(M,ns);
C = zeros(M,ns);
for i=1:M
    x_center = x(i);
    y_center = y(i);
    rd = sqrt((x_center-x).^2+(y_center-y).^2);
    [rd,ix] = sort(rd);
    F(i,:) = rd(1:ns);
    C(i,:) = ix(1:ns);
    ci(i) = rd(ns)*sqrt(ns)*rd(2)/(0.02*(sum(F(i,1:ns))));
end

%%
c = 4;
minK = 1e8;
maxK = 1e12;
dc = 0.01;

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
Re=800;
U=zeros(N,1);
V=zeros(N,1);
U(north)=1;
%% Pressure Equation
Lpsi=Dxx+Dyy;
Lpsi(north,:)=Dy(north,:);
Lpsi(south,:)=-Dy(south,:);
Lpsi(east,:)=Dx(east,:);
Lpsi(west,:)=-Dx(west,:);
%%
k=0;
I = eye(size(Dx));
for i=1:Nt-1
    k=k+1;
    
    %%
    RHS_u = (1/dt)*U ;
    RHS_v = (1/dt)*V ; 
    %%
    DDu = (1/dt).*I + U.*Dx + V.*Dy + (1/(Re))*(Dxx+Dyy) ;
    DDv = (1/dt).*I + U.*Dx + V.*Dy + (1/(Re))*(Dxx+Dyy) ;    
    %%
    U_star = pinv(DDu)*RHS_u ;
    V_star = pinv(DDv)*RHS_v ;
    
    %%
    RHS_Poisson = (1/dt)*(Dx*U_star + Dy*V_star);
    
    RHS_Poisson(Boundary) = 0;
    
    P = Lpsi\RHS_Poisson;
    %%
    U = U_star - dt*Dx*P;
    
    V = V_star - dt*Dy*P;
    
    U(Boundary)=0;  U(north)=1;
    
    V(Boundary)=0;
    
    if (k==TSCREEN)
        plot(0)
        hh=streamslice(X,Y,reshape(U,size(X)),reshape(V,size(X)),25,'k','noarrow');
        set( hh, 'Color', 'b');
        hold off
        drawnow
        k=0;
    end

    
    T(i+1)
end

