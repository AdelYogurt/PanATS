function [direction,radius,base_point,fit_error]=fitCylinder(point_list)
% newton iteration to fit the cylinder
%
point_number=size(point_list,1);

iter_max=10;
iter=0;
torl=1e-3;
base_point=mean(point_list,1); % a, b, c
radius=sum(sqrt(sum((point_list-base_point).^2,2)))/point_number; % r
direction=base_point;

warning off
while iter < iter_max
    P_PB=(point_list-base_point);
    P=sum(direction.*P_PB,2);
    C0=sum(P_PB.^2,2)-radius*radius-P.*P;
    dF_dPB=-2*P_PB+2*direction.*P;
    dF_dr=-2*ones(point_number,1)*radius;
    dF_dD=-2*P_PB.*P;

    Grad=[dF_dPB,dF_dr,dF_dD];
    Fval=C0;
    dVar=-Grad\Fval;
    base_point=base_point+dVar(1:3)';
    radius=radius+dVar(4);
    direction=direction+dVar(5:7)';
    direction=direction/norm(direction);

    if norm(dVar) < torl
        break;
    end

    iter=iter+1;
end
warning on

% calculate error
% fit_error=(sum((point_list-base_point).^2,2)-radius*radius-sum((direction.*(point_list-base_point)).^2,2))/(radius*radius);
% fit_error=sum(abs(fit_error))/point_number;

radius_sq=radius*radius;
radius_sq_real=sum((point_list-base_point).^2-(direction.*(point_list-base_point)).^2,2);
fit_error=mean(abs(radius_sq_real-radius_sq));
end

% function [direction,radius,base_point,fit_error]=fitCylinder(point_list)
% % newton iteration to fit the cylinder
% %
% point_number=size(point_list,1);
% 
% iter_max=10;
% iter=0;
% torl=1e-3;
% base_point=mean(point_list,1); % a, b, c
% radius=sum(sqrt(sum((point_list-base_point).^2,2)))/point_number; % r
% direction=base_point;
% 
% while iter < iter_max
%     % updata direction first
%     P_PB=(point_list-base_point);
%     P=sum(direction.*P_PB,2);
%     dF_dD=-2*P_PB.*P;
%     C0=sum(P_PB.^2,2)-radius*radius-P.*P;
%     ddirection=-dF_dD\C0;
%     direction=direction+ddirection';
%     direction=direction/norm(direction);
% 
%     % updata base_point and radius
%     P_PB=(point_list-base_point);
%     P=sum(direction.*P_PB,2);
%     dF_dPB=-2*P_PB+2*direction.*P;
%     dF_dr=-2*ones(point_number,1)*radius;
%     C0=sum(P_PB.^2,2)-radius*radius-P.*P;
%     dVar=-[dF_dPB,dF_dr]\C0;
%     base_point=base_point+dVar(1:3)';
%     radius=radius+dVar(4);
% 
%     if norm(dVar) < torl
%         break;
%     end
% 
%     iter=iter+1;
% end
% 
% % calculate error
% fit_error=(sum((point_list-base_point).^2,2)-radius*radius)-...
%     sum((direction.*(point_list-base_point)).^2,2)/(radius*radius);
% fit_error=sum(abs(fit_error))/point_number;
% end

% function [direction,radius,base_point,fit_error]=fitCylinder(point_list)
% % least squares to fit the cylinder
% %
% point_number=size(point_list,1);
% 
% normal=pca(point_list);
% direction=normal(:,1)';
% C=sum(direction.^2);
% 
% X=point_list(:,1);
% Y=point_list(:,2);
% Z=point_list(:,3);
% l=direction(1);
% m=direction(2);
% n=direction(3);
% 
% % notice l, m, n is known
% matrix=[-2*point_list,ones(point_number,1)];
% U=-sum(point_list.^2,2)+sum(l^2*X.^2 +2*l*m*X.*Y+2*l*n*X.*Z+m^2*Y.^2+2*m*n*Y.*Z +n^2*Z.^2,2)/C;
% coeff=(matrix\U); % x0, y0, z0, D
% 
% base_point=coeff(1:3)';
% x0=base_point(1);
% y0=base_point(2);
% z0=base_point(3);
% 
% D=sum(base_point.^2)+sum(-2*l^2*X*x0+l^2*x0^2-2*l*m*x0*Y+2*l*m*x0*y0-2*l*n*x0*Z+2*l*n*x0*z0-2*l*m*X*y0-2*m^2*Y*y0+m^2*y0^2-2*m*n*y0*Z+2*m*n*y0*z0-2*l*n*X*z0-2*m*n*Y*z0-2*n^2*Z*z0+n^2*z0^2)/C-coeff(4);
% D=abs(D); % here have some problem, if not all point are on the cylinder
% radius=sqrt(D);
% 
% % calculate error
% fit_error=(sum((point_list-base_point).^2,2)-radius*radius)-...
%     sum((direction.*(point_list-base_point)).^2,2)/C/(radius*radius);
% fit_error=sum(abs(fit_error))/point_number;
% end
