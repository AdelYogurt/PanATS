clc;
% clear all;
close all hidden;

addpath([pwd,'\src']);

global g_geometry g_Point g_Element g_Marker...
    ADtree_marker_element HATS_element_list...
    streamline_output inviscid_output heat_output viscid_output FEM_output...
    post_output

% displayModel('P',FREESTREAM_PRESSURE)
% displayMarker('WAVERIDER')

PostModel()

% [draw_X,draw_data]=getSYMdata(post_output.log_P_point_list);
% figure();
% plot(draw_X,log(draw_data/FREESTREAM_PRESSURE)/log(10)+1,'LineStyle','none','Marker','.')
% plot(draw_X,draw_data,'LineStyle','none','Marker','.')

load('blunt_data.mat')
[draw_X,draw_data]=getSYMdata(post_output.Q_point_list);
[draw_X,index_list]=sort(draw_X);
draw_X=draw_X-draw_X(1);
draw_data=draw_data(index_list);

figure();
sgtitle('Blunt Heat Verification');
line(draw_X/568.7e-3,draw_data/Q_0,'Marker','.');
line(x,Q_blunt,'LineStyle','none','Marker','o','Color','r');
legend('Estimate','EXP')
set(gca,'Xlim',[-0.1,1]);
set(gca,'Ylim',[-0.3,2.3]);
xlabel('x/L');
ylabel('q/q_0');

error_list=zeros(1,length(x));
for x_index=1:length(x)
    data=interpolationLinear(x(x_index),draw_X/568.7e-3,draw_data/Q_0);
    error_list(x_index)=abs(data-Q_blunt(x_index))/Q_blunt(x_index);
end

function data=interpolationLinear(x,x_list,data_list)
index=1;
while x > x_list(index+1)
    index=index+1;
end
d_data=(data_list(index+1)-data_list(index))/(x_list(index+1)-x_list(index));
data=data_list(index)+d_data*(x-x_list(index));
end

function PostModel()
% post model, add value to all point
%
global g_geometry g_Point g_Element g_Marker...
    ADtree_marker_element HATS_element_list...
    streamline_output inviscid_output heat_output viscid_output...
    post_output

marker_element_number=size(HATS_element_list,1);
dimension=3;

marker_point_index_list=[HATS_element_list.point_index_list];
[marker_point_index_list,~,mapping_list]=unique(marker_point_index_list);

% load data
inflow_direction_list=streamline_output.inflow_direction_list;
streamline_length_list=streamline_output.streamline_length_list;
Cp_list=inviscid_output.Cp_list;
P_list=inviscid_output.P_list;
dFn_list=inviscid_output.dFn_list;
dMn_list=inviscid_output.dMn_list;
Q_list=heat_output.Q_list;
Re_x_list=heat_output.Re_x_list;
dFs_list=viscid_output.dFs_list;
dMs_list=viscid_output.dMs_list;

Cp_point_list=zeros(length(marker_point_index_list),1);
P_point_list=zeros(length(marker_point_index_list),1);
dFn_point_list=zeros(length(marker_point_index_list),1);
dMn_point_list=zeros(length(marker_point_index_list),1);
Q_point_list=zeros(length(marker_point_index_list),1);
Re_x_point_list=zeros(length(marker_point_index_list),1);
dFs_point_list=zeros(length(marker_point_index_list),1);
dMs_point_list=zeros(length(marker_point_index_list),1);

add_time_list=zeros(length(marker_point_index_list),1);
for element_index=1:marker_element_number
    point_index_list=mapping_list(element_index*3-2:element_index*3);
    
    Cp_point_list(point_index_list,:)=Cp_point_list(point_index_list,:)+Cp_list(element_index);
    P_point_list(point_index_list,:)=P_point_list(point_index_list,:)+P_list(element_index);
    dFn_point_list(point_index_list,:)=dFn_point_list(point_index_list,:)+dFn_list(element_index);
    dMn_point_list(point_index_list,:)=dMn_point_list(point_index_list,:)+dMn_list(element_index);
    Q_point_list(point_index_list,:)=Q_point_list(point_index_list,:)+Q_list(element_index);
    Re_x_point_list(point_index_list,:)=Re_x_point_list(point_index_list,:)+Re_x_list(element_index);
    dFs_point_list(point_index_list,:)=dFs_point_list(point_index_list,:)+dFs_list(element_index);
    dMs_point_list(point_index_list,:)=dMs_point_list(point_index_list,:)+dMs_list(element_index);
    
    add_time_list(point_index_list,:)=add_time_list(point_index_list,:)+1;
end
Cp_point_list=Cp_point_list./add_time_list;
P_point_list=P_point_list./add_time_list;
dFn_point_list=dFn_point_list./add_time_list;
dMn_point_list=dMn_point_list./add_time_list;
Q_point_list=Q_point_list./add_time_list;
Re_x_point_list=Re_x_point_list./add_time_list;
dFs_point_list=dFs_point_list./add_time_list;
dMs_point_list=dMs_point_list./add_time_list;

P_1=g_geometry.P_1;

post_output.Cp_point_list=Cp_point_list;
post_output.P_point_list=P_point_list;
post_output.log_P_point_list=log(P_point_list/P_1)/log(10)+1;
post_output.dFn_point_list=dFn_point_list;
post_output.dMn_point_list=dMn_point_list;
post_output.Q_point_list=Q_point_list;
post_output.Re_x_point_list=Re_x_point_list;
post_output.dFs_point_list=dFs_point_list;
post_output.dMs_point_list=dMs_point_list;
post_output.marker_point_index_list=marker_point_index_list;
end

function [draw_X,draw_data]=getSYMdata(data_list)
% read data on symmetry
%
global g_Point post_output

marker_point_index_list=post_output.marker_point_index_list;
draw_X=[];
draw_data=[];
for point_index=1:length(marker_point_index_list)
   if g_Point(marker_point_index_list(point_index),2) < 1e-6
       draw_X=[draw_X;g_Point(marker_point_index_list(point_index),1)];
       draw_data=[draw_data;data_list(point_index)];
   end
end
end

function writeNode(file_name,permission,point_index_list)
% function to print point coordinate to file
% include point index(global) point coordinate
%
global g_Point
if isempty(permission)
    permission='w+';
end

if ~strcmp(file_name(end-3:end),'.dat')
    file_name=[file_name,'.dat'];
end

dimension=size(g_Point,2);

file_point=fopen(file_name,permission);

% fprintf g_Point of point_index_list
for point_index=1:size(point_index_list,1)
    point_unit=point_index_list(point_index);
    
    % point index(global)
    fprintf(file_point,'%d',point_unit);
    
    % point coordinate
    for dimension_index=1:dimension
        fprintf(file_point,' % 9.8e',g_Point(point_unit+1,dimension_index));
    end
    
    fprintf(file_point,'\n');
end
fclose(file_point);
clear('file_point');
end
function writeMarker(file_name,permission,marker)
% write marker to file
% only support su2 format
% marker is cell include marker name, element number, g_Element
% write file is dat format
%
if isempty(permission)
    permission='w+';
end
if ~strcmp(file_name(end-3:end),'.dat')
    file_name=[file_name,'.dat'];
end

file_marker=fopen(file_name,permission);

% fprintf marker name
fprintf(file_marker,'MARKER_TAG= %s\n',marker{1});

% fprintf element number
fprintf(file_marker,'MARKER_ELEMS= %d\n',marker{2});

% fprintf g_Element
g_Element=marker{3};
for element_index=1:size(g_Element,1)
    element=g_Element(element_index,:);
    element_type=element(1);
    
    % element type
    fprintf(file_marker,'%d',element_type);
    
    % element point index
    switch(element_type)
        case 3
            point_number=2;
        case 5
            point_number=3;
        case 9
            point_number=4;
        case 10
            point_number=4;
        case 12
            point_number=8;
        case 13
            point_number=6;
        case 14
            point_number=5;
    end
    for point_index=1:point_number
        fprintf(file_marker,' %d',element(1+point_index));
    end
    fprintf(file_marker,'\n');
end
fclose(file_marker);
clear('file_airfoil');
end

function point_index_list=getElementNode(g_Element)
% get all point of g_Element
% avoid repeat and sort by global index
%
point_index_list=[];

for element_index=1:size(g_Element,1)
    element=g_Element(element_index,:);
    element_type=element(1);
    
    % element point index
    switch(element_type)
        case 3
            point_number=2;
        case 5
            point_number=3;
        case 9
            point_number=4;
        case 10
            point_number=4;
        case 12
            point_number=8;
        case 13
            point_number=6;
        case 14
            point_number=5;
    end
    point_index_list=[point_index_list;element(2:1+point_number)'];
end

% delete repeat and sort
point_index_list=unique(point_index_list);
end
function Marker_Element=getMarkerElement(marker_name)
% return specified marker
% marker is include all element unit of marker
% element include element type and point index
%
global g_Marker
Marker_Element=[];
for marker_index=1:size(g_Marker,1)
    if strcmp(g_Marker{marker_index,1},marker_name)
        Marker_Element=g_Marker{marker_index,3};
    end
end
if isempty(Marker_Element)
    error('getMarkerElement: no marker found')
end
end
