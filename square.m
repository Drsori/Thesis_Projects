

x = -1:0.2:1;
y = -1:0.2:1;

[Y,X] = meshgrid(x,y);
B = find(X==-1|X==1|Y==-1|Y==1);

plot(X(B),Y(B),'ro')