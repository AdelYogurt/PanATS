clc;
clear;
close all hidden;

% xc = 1.21;
% yc = 2.32;
% zc = 4.32;
% 
% xr = 2.78;
% yr = 5.76;
% zr = 1.51;
% 
% [x,y,z] = ellipsoid(xc,yc,zc,xr,yr,zr,20);
% 
% plot3(x(:),y(:),z(:),'.');
% 
% x=x(:);
% y=y(:);
% z=z(:);
% 
% X = [x.*x y.*y z.*z x.*y x.*z y.*z x y z];
% Y = ones(length(x),1);
% 
% C = inv(X'*X)*X'*Y;
% 
% M=[C(1) C(4)/2 C(5)/2;
%     C(4)/2 C(2) C(6)/2;
%     C(5)/2 C(6)/2 C(3)];
% 
% Cent = -0.5*[C(7) C(8) C(9)]*inv(M)
% 
% S = Cent*M*Cent'+1;
% [U,V]=eig(M);
% 
% [~,indx] = max(abs(U(1,:)));
% [~,indy] = max(abs(U(2,:)));
% [~,indz] = max(abs(U(3,:)));
% 
% R = [sqrt(S/V(indx,indx)) sqrt(S/V(indy,indy)) sqrt(S/V(indz,indz))]


point=[2,0,1];
point_ref=[0,2,1];
point_start=[0,0,0];
surface_flow=[0.5,0.5,0.5];

point_end=calCrossPoint(point,point_ref,point_start,surface_flow)

function point_end=calCrossPoint(point,point_ref,point_start,surface_flow)
% calculate cross point coordinate as surface_flow downstream point_end
%
surface_flow=surface_flow/norm(surface_flow);
dr_edge=point_ref-point;
dr_edge_length=norm(dr_edge);
dr_edge_norm=dr_edge/dr_edge_length;
point_proj_length=norm(cross(surface_flow,point-point_start));
point_end=point+point_proj_length/norm(cross(surface_flow,dr_edge_norm))/norm(dr_edge)*dr_edge;
end

