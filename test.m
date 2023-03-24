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

% judgeOriginSurround([0,0;1,0;1,1;0,1],1e-6)

load('geo.mat')
judgeOriginSurround(point_2D_list(1:end,[1,2]),0.01)

function surround_flag=judgeOriginSurround(curve_origin,geometry_torlance)
% function to judge if origin point(0,0) is surrounded by line
% calculate curve cross positive axis x times, if is odd number, surround
%
surround_flag=0;

% shifting line with geometry_torlance
curve=curve_origin;
line(curve_origin(:,1),curve_origin(:,2))
for edge_index=1:size(curve_origin,1)
    point_1=curve_origin(edge_index,:);
    if edge_index == size(curve_origin,1)
        point_2=curve_origin(1,:);
        edge_index_next=1;
    else
        point_2=curve_origin(edge_index+1,:);
        edge_index_next=edge_index+1;
    end

    A=point_1(2)-point_2(2);
    B=point_2(1)-point_1(1);
    move_dir=-[A,B]/norm([A,B]);

    % move point of edge
    curve(edge_index,:)=curve(edge_index,:)+move_dir*geometry_torlance;
    curve(edge_index_next,:)=curve(edge_index_next,:)+move_dir*geometry_torlance;
    line(curve(:,1),curve(:,2))
end

cross_time=0;
for edge_index=1:size(curve,1)
    point_1=curve(edge_index,:);
    if edge_index == size(curve,1)
        point_2=curve(1,:);
    else
        point_2=curve(edge_index+1,:);
    end

    if (((point_1(2) < -0) && (point_2(2) < -0)) || ...
            ((point_1(2) > 0) && (point_2(2) > 0)))
        % if edge do not cross axis x
        continue;
    end

    % edge cross axis x, calculate cross point
    A=point_1(2)-point_2(2);
    C=point_1(1)*point_2(2)-point_2(1)*point_1(2);
    cross_point=-C/A;

    % check if edge is horizontal
    if abs(A) > 0 % no horizontal
        if cross_point > -0
            cross_time=cross_time+1;
        end
    else % horizontal
        if ((point_1(1) > -0) ||...
                (point_2(1) > -0))
            % cross point must in line range
            % line overlap X axis
            cross_time=cross_time+1;
        end
    end
end
if (mode(cross_time,2) == 1)
    surround_flag=1;
end
end
