





clear 
clc
format shorte
close all
tic
%%
a = 0;  b = 1;
np = 40;
c = 0.5;

h = (b-a)/np;
x = a:h:b;
y = a:h:b;
[X,Y] = meshgrid(x);
X = X(:);  Y = Y(:);
N = length(X);

%% Without corner
Boundary = find( X==a | X==b | Y==a | Y==b ...
    & not (X==a & Y==b) ...
    & not (X==a & Y==a)...
    & not (X==b & Y==a)...
    & not (X==b & Y==b));
East = find( X==a & Y~=a & Y~=b );
West = find( X==b & Y~=a & Y~=b );
South = find( Y==a  & X~=a & X~=b);
North = find( Y==b & X~=a & X~=b );

%% With corner

% Boundary = find( X==a | X==b | Y==a | Y==b);
% 
% East = find( X==a );
% West = find( X==b );
% South = find( Y==a );
% North = find( Y==b );


%%
T = 1;
tau = 1e-1;
Nt = T/tau;
t = 0:tau:T;
%%
u = @(x,y,t) exp(-t).*sin(pi*x).*sin(pi*y);
v = @(x,y,t) exp(-2*t).*sin(2*pi*x).*sin(2*pi*y);

f = @(x,y,t) -exp(-t).*sin(pi*x).*sin(pi*y) + 2*(pi^2)*exp(-t).*sin(pi*x).*sin(pi*y)-exp(-2*t).*sin(2*pi*x).*sin(2*pi*y);
% g = @(x,y,t) - 2*exp(-t).*sin(pi*x).*sin(pi*y) - 2*exp(-2*t).*sin(2*pi*x).*sin(2*pi*y);
g = @(x,y,t) -2*exp(-t).*sin(pi*x).*sin(pi*y) + (8*pi^2-2).*exp(-2.*t).*sin(2*pi.*x).*sin(2*pi.*y);

%% Derivatives of exact solution
    ux=@(x,y,t) pi*exp(-t).*cos(pi.*x).*sin(pi*y);
    uy=@(x,y,t) pi*exp(-t).*cos(pi.*y).*sin(pi*x);
    
    vx=@(x,y,t) 2*pi*exp(-2*t).*cos(2*pi.*x).*sin(2*pi*y);
    vy=@(x,y,t) 2*pi*exp(-2*t).*cos(2*pi.*y).*sin(2*pi*x);

%% RBFs

phi    = @(r,c) exp(-c^2.*r.^2);
dxphi  = @(r,rx,c) -2*c^2*rx.*exp(-c^2.*r.^2);
dyphi  = @(r,ry,c) -2*c^2*ry.*exp(-c^2.*r.^2);
dxxphi = @(r,rx,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(rx.^2));
dyyphi = @(r,ry,c) 2*c^2*exp(-(c*r).^2).*( -1 + 2*(c^2).*(ry.^2));


rx = X-X';
ry = Y-Y';
r = sqrt( (rx).^2 + (ry).^2 );
%% 
A = phi(r,c);

Ax = dxphi(r,rx,c);
Ay = dyphi(r,ry,c);

Axx = dxxphi(r,rx,c);
Ayy = dyyphi(r,ry,c);

D = Axx + Ayy;
%%
I = eye(N);
Z =zeros(N);

A1 = A -0.5*tau*D;
A2 = -0.5*tau*A;
A3 = -tau*A ;
A4 = A - 0.5*tau*D;

%%
A1(East,:)  = -Ax(East,:);
A1(West,:)  =  Ax(West,:);
A1(South,:) = -Ay(South,:);
A1(North,:) =  Ay(North,:);

A4(East,:)  = -Ax(East,:);
A4(West,:)  =  Ax(West,:);
A4(South,:) = -Ay(South,:);
A4(North,:) =  Ay(North,:);

%   A2(Boundary,:) = I(Boundary,:);
%   A3(Boundary,:) = I(Boundary,:);
%%
H_left = [A1, A2 ; A3, A4];


H_right = [A + 0.5*tau*D, 0.5*tau*A;...
            tau*A,  0.5*tau*D];

H=pinv(H_left);
%% Initial value
   U = pinv(A)*u(X,Y,0);
   V = pinv(A)*v(X,Y,0);
    
UV = [U ; V];

for n=1:Nt
    time = (n-0.5)*tau;
    F = H_right*UV + tau*[f(X,Y,time);g(X,Y,time)];
    
    F(East,:)  =  -ux(X(East),Y(East),time);
    F(West,:)  =   ux(X(West),Y(West),time);
    F(South,:) =  -uy(X(South),Y(South),time);
    F(North,:) =   uy(X(North),Y(North),time);
   
    
    F(East+N,:)  =  -vx(X(East),Y(East),time);
    F(West+N,:)  =   vx(X(West),Y(West),time);
    F(South+N,:) =  -vy(X(South),Y(South),time);
    F(North+N,:) =   vy(X(North),Y(North),time);
    
    UV = H*F;
end

U = A*UV(1:N);
V = A*UV(N+1:end);

Uex=u(X,Y,T);
Vex=v(X,Y,T);

Error_inf_U=norm(U-Uex,inf)
Error_inf_V=norm(V-Vex,inf)





