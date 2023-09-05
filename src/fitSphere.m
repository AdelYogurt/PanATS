function [center_point,radius,fit_error]=fitSphere(point_list)
% least squares to fit the sphere
%
point_number=size(point_list,1);

matrix=[-2*point_list,ones(point_number,1)];
U=-sum(point_list.^2,2);
coeff=(matrix\U); % xc, yc, zc, D
center_point=coeff(1:3);
% D=xc*xc+yc*yc+zc*zc-r*r
radius=sqrt(sum(center_point.^2)-coeff(4));

% calculate error
% fit_error=(sum((point_list-center_point').^2,2)-radius*radius)/(radius*radius);
% fit_error=sum(abs(fit_error))/point_number;

radius_sq=radius*radius;
radius_sq_real=sum((point_list-center_point').^2,2);
fit_error=mean(abs(radius_sq_real-radius_sq));
end