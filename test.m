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

curve=[0,0;2,0;1,2;-2,1];
line(curve(:,1),curve(:,2));
curve=geoCurveOffset([0,0;2,0;1,2;-2,1],0.5);
line(curve(:,1),curve(:,2));
axis equal;

function curve=geoCurveOffset(curve,offset)
% shifting line with offset
% direction is outside
%
edge_number=size(curve,1);

% calculate line data and move line
edge_data_list=zeros(edge_number,3); % A, B, C
for edge_index=1:edge_number
    point_1=curve(edge_index,:);
    if edge_index == size(curve,1)
        point_2=curve(1,:);
    else
        point_2=curve(edge_index+1,:);
    end

    dr=point_2-point_1;
    
    edge_data_list(edge_index,1)=-dr(2); % A (dr_y=-A)
    edge_data_list(edge_index,2)=dr(1); % B (dr_x=B)

    % move line
    edge_data_list(edge_index,3)=point_1(1)*point_2(2)-point_1(2)*point_2(1)+...
        offset/norm(dr)*sum(dr.^2); % C
end

% solve new cross point
for edge_index=1:edge_number
    if edge_index == 1
        edge_prev_index=edge_number;
    else
        edge_prev_index=edge_index-1;
    end

    matrix=[edge_data_list(edge_index,1),edge_data_list(edge_index,2);
        edge_data_list(edge_prev_index,1),edge_data_list(edge_prev_index,2)];
    curve(edge_index,:)=matrix\[-edge_data_list(edge_index,3);-edge_data_list(edge_prev_index,3)];
end

end
