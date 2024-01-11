function [center_point,radius]=fitEllipsoid(point_list)
% least squares to fit the cylinder ellipsoid
%
x=point_list(:,1);
y=point_list(:,2);
z=point_list(:,3);

D=[x.*x,y.*y,z.*z,x.*y,x.*z,y.*z,x,y,z];
P=ones(length(x),1);

C=D\P;

M=[C(1) C(4)/2 C(5)/2;
    C(4)/2 C(2) C(6)/2;
    C(5)/2 C(6)/2 C(3)];

center_point=-0.5*M\[C(7);C(8);C(9)];

S=center_point'*M*center_point+1;
lamada=eig(M);

radius=[sqrt(S/lamada(1));sqrt(S/lamada(2));sqrt(S/lamada(3))];

end