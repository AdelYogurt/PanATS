clc;
clear;
close all hidden;

jacobian=[
    6,0;0,2];
bias=[0;0];

for j=-0.5:0.2:0.5
    for i=-0.5:0.1:0.5
        [t,y]=ode45(@(t,y) nolinear(t,y,jacobian,bias),[0,0.1],[i,j]);
        line(y(:,1),y(:,2));
    end
end

point_list=[
    -0.5,-0.5,0;
    0.5,-0.5,0;
    0.5,0.5,0;
    -0.5,0.5,0;];

surface_flow_list=point_list(:,1:2)*jacobian'+bias';
surface_flow_list=[surface_flow_list,zeros(size(surface_flow_list,1),1)]

function dy=nolinear(t,y,jacobian,bias)
dy=jacobian*y+bias;
end
