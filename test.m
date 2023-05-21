clc;
clear;
close all hidden;

% syms X x0 Y y0 Z z0 l m n C
% fcn=(l*(X-x0)+m*(Y-y0)+n*(Z-z0))^2;
% expand(fcn)

load('matlab.mat');

scatter3(point_arho_list(:,1),point_arho_list(:,2),point_arho_list(:,3));
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');

[direction,radius,base_point,error]=fitCylinder(point_arho_list)


% x = [0 1-1 0 0]';
% y = [0 0 0 1-1]';
% z = [1 0 0 0 0]';
% point_list = [x y z];
% [center_point,r]=fitSphere(point_list)

function [direction,radius,base_point,error]=fitCylinder(point_list)
% least squares to solve the cylinder
%
point_number=size(point_list,1);

normal=pca(point_list);
direction=normal(:,1)';
C=sum(direction.^2);

X=point_list(:,1);
Y=point_list(:,2);
Z=point_list(:,3);
l=direction(1);
m=direction(2);
n=direction(3);

% notice l, m, n is known
matrix=[-2*point_list,ones(point_number,1)];
U=-sum(point_list.^2,2)+sum(l^2*X.^2 +2*l*m*X.*Y+2*l*n*X.*Z+m^2*Y.^2+2*m*n*Y.*Z +n^2*Z.^2,2)/C;
coeff=(matrix\U); % x0, y0, z0, D

base_point=coeff(1:3)';
x0=base_point(1);
y0=base_point(2);
z0=base_point(3);

radius=sqrt(sum(base_point.^2)+sum(-2*l^2*X*x0+l^2*x0^2-2*l*m*x0*Y+2*l*m*x0*y0-2*l*n*x0*Z+2*l*n*x0*z0-2*l*m*X*y0-2*m^2*Y*y0+m^2*y0^2-2*m*n*y0*Z+2*m*n*y0*z0-2*l*n*X*z0-2*m*n*Y*z0-2*n^2*Z*z0+n^2*z0^2)/C-coeff(4));
radius=real(radius);

% calculate error
error_list=(sum((point_list-base_point).^2,2)-radius*radius)-...
    sum((direction.*(point_list-base_point)).^2,2)/C;
error=sum(error_list.^2)/point_number;
end

