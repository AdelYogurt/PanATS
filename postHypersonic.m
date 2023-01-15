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

function [T,P,rou,a,g]=getAtmosphereEnvironment(Z)
% base on altitude calculate atmosphere parameter
% Z is altitude m
% return temperature pressure density speed_of_sound acceleration_of_gravity
%
Z=Z/1e3;

R_earth=    6.356766e3;     % earth ratio km
T_SL=       2.8815e2;       % sea level temperature
P_SL=       1.01325e5;      % sea level pressure
rou_SL=     1.2250;         % sea level density
g_SL=       9.80665;        % sea level acceleration of gravity
a_SL=       3.40294e2;      % sea level speed of sound

% geopotential height H
H=Z/(1+Z/R_earth);

% calculate parameter
if Z <= 11.0191
    W=1-H/44.3308;
    T=W*T_SL;
    P=W^5.2553*P_SL;
    rou=W^4.2559*rou_SL;
    
elseif Z <= 20.0631
    W=exp((14.9647-H)/6.3416);
    T=216.650;
    P=1.1953*1e-1*W*P_SL;
    rou=1.5898*1e-1*W*rou_SL;
    
elseif Z <= 32.1619
    W=1+(H-24.9021)/221.552;
    T=221.552*W;
    P=2.5158*1e-2*W^-34.1629*P_SL;
    rou=3.2722*1e-2*W^-35.1829*rou_SL;
    
elseif Z <= 47.3501
    W=1+(H-39.4799)/89.4107;
    T=250.350*W;
    P=2.8338*1e-3*W^-12.2011*P_SL;
    rou=3.2617*1e-3*W^-13.2011*rou_SL;
    
elseif Z <= 51.4125
    W=exp((48.6252-H)/7.9223);
    T=270.650;
    P=8.9155*1e-4*W*P_SL;
    rou=9.4920*1e-4*W*rou_SL;
    
elseif Z <=  71.8020
    W=1-(H-59.4390)/88.2218;
    T=247.021*W;
    P=2.1671*1e-4*W^12.2011*P_SL;
    rou=2.5280*1e-4*W^11.2011*rou_SL;
    
elseif Z <= 86.0000
    W=1-(H-78.0303)/100.2950;
    T=200.590*W;
    P=1.2274*1e-5*W^17.0816*P_SL;
    rou=1.7632*1e-5*W^16.0816*rou_SL;
    
elseif Z <= 91.0000
    W=exp((87.2848-H)/5.4700);
    T=186.870;
    P=(2.2730+1.042*1e-3*H)*1e-6*W*P_SL;
    rou=3.6411*1e-6*W*rou_SL;
    
elseif Z <= 150.0000
    T=-0.0040*Z^3+1.5054*Z^2-177.5620*Z+6.8929*1e-3;
    P=2.532*1e6*exp(-0.1829*Z)+0.1403*exp(-0.03698*Z);
    rou=70.22*1e6*exp(-0.1874*Z)+1.734*1e-5*exp(-0.05828*Z);
    
else
    
    
end

a=20.0468*sqrt(T);
g=g_SL/(1+Z/R_earth)^2;

end

function displayMarker(marker_name)
% draw marker point
%
global g_Point

dimension=size(g_Point,2);

marker_g_Element=getMarkerElement(marker_name);

if dimension == 2
    %     line(g_Point(point_index_list+1,1),g_Point(point_index_list+1,2),...
    %         'LineStyle','none','Marker','o');
    for element_index=1:size(marker_g_Element,1)
        element=marker_g_Element(element_index,:);
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
        
        element_point_index_list=element(2:1+point_number);
        line(g_Point(element_point_index_list+1,1),g_Point(element_point_index_list+1,2));
        
    end
    xlabel('x');
    ylabel('y');
elseif dimension == 3
    %     line(g_Point(point_index_list+1,1),g_Point(point_index_list+1,2),g_Point(point_index_list+1,3),...
    %         'LineStyle','none','Marker','o');
    for element_index=1:size(marker_g_Element,1)
        element=marker_g_Element(element_index,:);
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
        
        element_point_index_list=element(2:1+point_number);
        line(g_Point(element_point_index_list+1,1),g_Point(element_point_index_list+1,2),g_Point(element_point_index_list+1,3));
        
    end
    xlabel('x');
    ylabel('y');
    zlabel('z');
end
axis equal;

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
