clc;
clear;
close all hidden;

a1=1;a2=0;
b1=-5;b2=1;
c1=-2;c2=3;

xcp=(a2*c1-a1*c2)/(b1*c2-b2*c1)
ycp=(a1*b2-a2*b1)/(b1*c2-b2*c1)

[eig_vector,eig_value]=eig([b1,c1;b2,c2])
lamada=eig_value(1,1);
miu=eig_value(2,2);

equation=['Dx=',num2str(a1),'+',num2str(b1),'*x',num2str(c1),'*y',',',...
    'Dy=',num2str(a2),'+',num2str(b2),'*x',num2str(c2),'*y'];

syms x(t) y(t)
equation = [diff(x,t)==a1+b1*x+c1*y, diff(y,t)==a2+b2*x+c2*y];

% draw
for x_initial=-2:1:2
    cond = [x(0)==x_initial,y(0)==1];
    [xSol(t),ySol(t)] = dsolve(equation,cond);
    line(xSol(-1:0.1:1),ySol(-1:0.1:1));
end

line(xcp,ycp,'Marker','o');

x_list=-5:0.5:5;
y_list=-5:0.5:5;
hold on;
for x_index=1:length(x_list)
    x=x_list(x_index);
    for y_index=1:length(y_list)
        y=y_list(y_index);
        quiver(x,y,a1+b1*x+c1*y,a2+b2*x+c2*y,'AutoScaleFactor',0.1);
    end
end
hold off

axis equal;

% load('matlab.mat');
% % point=[
% %     0,0;
% %     1,0;
% %     0,1;];
% 
% % velocity=[
% %     2,-0.5;
% %     1,-0.2;
% %     2,0.5;];
% 
% point=[0,0;d_center_point_list(:,1:2)];
% velocity=[abs_Ve,0;Ve_list(:,1:2)];
% 
% order=[1,2,3,4,1];
% hold on;
% line(point(order,1),point(order,2));
% quiver(point(1,1),point(1,2),velocity(1,1),velocity(1,2));
% quiver(point(2,1),point(2,2),velocity(2,1),velocity(2,2));
% quiver(point(3,1),point(3,2),velocity(3,1),velocity(3,2));
% quiver(point(4,1),point(4,2),velocity(4,1),velocity(4,2));
% hold off
% 
% matrix=reshape([d_center_point_list,zeros(element_nearby_number,6),d_center_point_list]',6,element_nearby_number*2)';
% matrix=[1,0,0,0,0,0;0,0,0,0,0,1;matrix]
% velocity=velocity';
% colume=Ve_list(:)
% 
% matrix=[
%     1,point(1,:),zeros(1,3);
%     zeros(1,3),1,point(1,:);
%     1,point(2,:),zeros(1,3);
%     zeros(1,3),1,point(2,:);
%     1,point(3,:),zeros(1,3);
%     zeros(1,3),1,point(3,:);
%     1,point(4,:),zeros(1,3);
%     zeros(1,3),1,point(4,:)]
% colume=[
%     velocity(1,:)';
%     velocity(2,:)';
%     velocity(3,:)';
%     velocity(4,:)']
% 
% coefficient=matrix\colume;
% 
% A=[coefficient(1);coefficient(4)]
% BC=[coefficient(2:3)';coefficient(5:6)']
% 
% -BC\A
% 
% [V,D]=eig(BC)
% 
% a1=coefficient(1);
% b1=coefficient(2);
% c1=coefficient(3);
% a2=coefficient(4);
% b2=coefficient(5);
% c2=coefficient(6);
% 
% equation=['Dx=',num2str(a1),'+',num2str(b1),'*x',num2str(c1),'*y',',',...
%     'Dy=',num2str(a2),'+',num2str(b2),'*x',num2str(c2),'*y'];
% 
% syms x(t) y(t)
% equation = [diff(x,t)==a1+b1*x+c1*y, diff(y,t)==a2+b2*x+c2*y];
% 
% % draw
% for x_initial=(-2:1:2)*1e-3
%     cond = [x(0)==x_initial,y(0)==0];
%     [xSol(t),ySol(t)] = dsolve(equation,cond);
%     line(xSol((-1:0.1:1)*1e-3),ySol((-1:0.1:1)*1e-3));
% end
% 
% x_list=(-5:0.5:5)*1e-3;
% y_list=(-5:0.5:5)*1e-3;
% hold on;
% for x_index=1:length(x_list)
%     x=x_list(x_index);
%     for y_index=1:length(y_list)
%         y=y_list(y_index);
%         quiver(x,y,a1+b1*x+c1*y,a2+b2*x+c2*y,'AutoScaleFactor',0.001);
%     end
% end
% hold off
% 
% axis equal;