
clear all
clc
close all
format shorte
%% Collocation Points
n = 20;
c = 8;
a = 0;
b = 1;
h = (b-a)/n;
x = a:h:b;
[X,Y] = meshgrid(x);
x = X(:);   y = Y(:);
N = length(x);
%%
N = length(x);
Boundary=find(x==a | x==b | y==a | y==b);
Interior=find(x~=a | x~=b | y~=a | y~=b);
%%
west=find(x==a & y~=a & y~=b);
east=find(x==b & y~=a & y~=b);
south=find(y==a & x~=a & x~=b);
north=find(y==b & x~=a & x~=b);
%%
% border1=[west;east;south;north];
% corner=(1:N)';
% corner([border1;Interior])=[];
%% Time discrete scheme
dt=1e-3;
T_final=10;
T=0:dt:T_final;
Nt=length(T);
TSCREEN=10;
%% Basis Functions and Required Matrices
Phi = @(r,c) exp(-(c*r).^2);
drPhi = @(r,rx,c) -2*c^2*rx.*exp(-(c*r).^2);
d2rPhi = @(r,rx,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rx.^2));

rx = x-x';
ry = y-y';

D = sqrt( rx.^2 + ry.^2 );
A = Phi(D,c);

Ax = drPhi(D,rx,c);
Ay = drPhi(D,ry,c);

Axx = d2rPhi(D,rx,c);
Ayy = d2rPhi(D,ry,c);
%%
I=eye(N);
Dx = Ax/A;
Dy = Ay/A;
Dxx = Axx/A;
Dyy = Ayy/A;
%%
Re=10;
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
for i=1:Nt-1
    k=k+1;
    
    %%
    RHS_u = (1/dt)*U ;
    RHS_v = (1/dt)*V ; 
    %%
    DDu = (1/dt).*I + U.*Dx + V.*Dy - (1/(Re))*(Dxx+Dyy) ;
    DDv = (1/dt).*I + U.*Dx + V.*Dy - (1/(Re))*(Dxx+Dyy) ;    
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
%        Go back in real space omega in real space for plotting
        plot(0)
        hh=streamslice(X,Y,reshape(U,size(X)),reshape(V,size(X)),25,'k','noarrow');
        set( hh, 'Color', 'b');
        hold off
        drawnow
        k=0;
    end
    T(i+1)
end

