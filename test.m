clc;
clear;
close all hidden;

xc = 1.21;
yc = 2.32;
zc = 4.32;

xr = 2.78;
yr = 5.76;
zr = 1.51;

[x,y,z] = ellipsoid(xc,yc,zc,xr,yr,zr,20);

plot3(x(:),y(:),z(:),'.');

x=x(:);
y=y(:);
z=z(:);

X = [x.*x y.*y z.*z x.*y x.*z y.*z x y z];
Y = ones(length(x),1);

C = inv(X'*X)*X'*Y;

M=[C(1) C(4)/2 C(5)/2;
    C(4)/2 C(2) C(6)/2;
    C(5)/2 C(6)/2 C(3)];

Cent = -0.5*[C(7) C(8) C(9)]*inv(M)

S = Cent*M*Cent'+1;
[U,V]=eig(M);

[~,indx] = max(abs(U(1,:)));
[~,indy] = max(abs(U(2,:)));
[~,indz] = max(abs(U(3,:)));

R = [sqrt(S/V(indx,indx)) sqrt(S/V(indy,indy)) sqrt(S/V(indz,indz))]