clc;
% clear;
close all hidden;

global user_model
geometry_torlance=1e-12;

user_model=preModelCFG('INP_plate.cfg');
preModelPanel();
preModelStreamline();
point_list=[
    0,0,0;
    1,0,0;
    1,1,0;
    0,1,0;];
marker.name='Part-plate';
marker.element_number=1;
marker.element_list=HATSElement(5,[1,2,3,4]);
marker_list=marker;

user_model.point_list=point_list;
user_model.marker_list=marker_list;

preModelPanel()
surface_flow_list=[
    -0.0    1.0         0
    0.0    1.0         0
    0.0135    0.99         0
    -0.0135    0.99         0];
user_model.surface_flow_list=surface_flow_list;
MARKER_MONITERING=user_model.MARKER_MONITORING;



% user_model=preModelCFG('INP_plate.cfg');
% preModelPanel();
% preModelStreamline();
% 
% point_list=user_model.point_list;
% marker_list=user_model.marker_list;
% 
% vertex_list=user_model.vertex_list;
% MARKER_MONITERING=user_model.MARKER_MONITORING;
% 
% % calculate surface velocity
% surface_flow_list=[zeros(size(point_list,1),4)]; % data, repeat times
% for point_index=1:size(point_list,1)
%     point=point_list(point_index,:);
%     surface_flow_list(point_index,1:3)=[point(1)-0.5,(point(2)-0.5),0];
% end
% 
% user_model.surface_flow_list=surface_flow_list;

% identift specify line
for marker_moniter_index=1:size(MARKER_MONITERING)
    marker_element_list=getMarkerElement(MARKER_MONITERING(marker_moniter_index),marker_list);
    for element_index=1:size(marker_element_list)
        element=marker_element_list(element_index);

        % load element point data
        point_index_list=element.point_index_list;
        point_center=sum(point_list(point_index_list,1:3))/length(point_index_list);
        Ve_list=surface_flow_list(point_index_list,1:3);
        point_2D_list=(point_list(point_index_list,1:3));
       
%         if element_index == 90
%             disp('');
%         end

        normal_vector=element.normal_vector;

        if judgeSpecify(normal_vector,point_2D_list,Ve_list)
            element.stagnation=true(1);
        end 
    end
end

displayMarker('Part-plate')


% load("matlab.mat");
% judgeSpecify(normal_vector,point_arou_list,Ve_list)


function cross_flag=judgeSpecify(normal_vector,point_list,Ve_list)
topology_torlance=1e-6;
precision_torlance=1e-12;

cross_flag=0;

% project to 2D
e3=normal_vector;
e1=point_list(1,1:3)-point_list(2,1:3);
e1=e1-(e1*e3')*e3;
e1=e1/norm(e1);
e2=cross(e3,e1);

% transform local cooridinate to 2D cooridinate
% all vector just simply project to element panel
% Ve_list=Ve_list*[e1;e2;e3]';
% point_2D_list=point_list*[e1;e2;e3]';
point_2D_list=point_list;
point_2D_list(:,end)=1;

% interpolation fit flow
coefficient=(point_2D_list\Ve_list(:,[1,2]))';

jacobian=coefficient(:,1:2); % [b1,c1;b2,c2]
jacobian(abs(jacobian) < precision_torlance)=0; % solve precision
bias=coefficient(:,3) % [a1;a2]

% jacobian matrix is negative, process stop
det_jacobian=(jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1));
if (((jacobian(1,1)+jacobian(2,2))^2-...
        4*det_jacobian) < -topology_torlance)
    return;
end

% evaluate eigenvalues of jacobian matrix
% if one of eigenvalues is zero, process stop
eig_value=eig(jacobian);
if sum(abs(eig_value) <= topology_torlance)
    return;
end

stagnation_point=-jacobian\bias;
% stagnation cannot to far 
if max(abs(stagnation_point)) > 10
    return;
end

% project element to canonical coordinates
[eig_vector,eig_value]=eig(jacobian);
eig_value=[eig_value(1,1),eig_value(2,2)];
if (eig_value(1) < eig_value(2))
    eig_value=fliplr(eig_value);
    eig_vector=fliplr(eig_vector);
end
% normalize eig_vector
point_2D_list=(point_2D_list(:,1:2)-stagnation_point')/eig_vector';
if det(eig_vector) < 0
    point_2D_list=flipud(point_2D_list);
end

% if (0,0) inside element
if judgeOriginSurround(point_2D_list(1:end,[1,2]),topology_torlance)
    cross_flag=1;
    return;
end

% judge phase portrait
% eigenvalue less than zero is concentrate
% X corresponds to the eigenvectors 1
% Y corresponds to the eigenvectors 2
if ((eig_value(1) > topology_torlance) && ...
        (eig_value(2) > topology_torlance))
    % repelling node, check small eigenvalue corresponded axis
    % if cross Y
    if judgeCrossY(point_2D_list(1:end,[1,2]),topology_torlance)
        cross_flag=1;
    end
elseif ((eig_value(1) < topology_torlance) && ...
        (eig_value(2) < topology_torlance))
    % attracting node, check
    % concentrate line is not stagnation point
    
elseif ((eig_value(1) > topology_torlance) && (eig_value(2) < -topology_torlance))
    % saddle, judge axis X and axis Y which is separation line
    % means which eigenvalue is large than zero
    % repelling node, check small eigenvalue corresponded axis
    if judgeCrossY(point_2D_list(1:end,[1,2]),topology_torlance)
        cross_flag=1;
    end
end

end

function setElementNoStagnation(element_list,element_number)
for element_index=1:element_number
    if isempty(element_list(element_index).stagnation)
        element_list(element_index).stagnation=false(1);
    elseif ~element_list(element_index).stagnation
        element_list(element_index).stagnation=false(1);
    end
end
end

function setElementStagnation(element_list,element_number)
for element_index=1:element_number
    element_list(element_index).stagnation=true(1);
end
end

function surround_flag=judgeOriginSurround(curve_origin,geometry_torlance)
% function to judge if origin point(0,0) is surrounded by line
% calculate curve cross positive axis x times, if is odd number, surround
%
surround_flag=0;

% shifting line with geometry_torlance
curve=curve_origin;
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

function cross_flag=judgeCrossX(curve,geometry_torlance)
% function to judge if curve cross X
%
cross_flag=0;
node_number=size(curve,1);
for node_index=1:node_number
    node_next_index=node_index+1;
    if node_next_index > node_number
        node_next_index=1;
    end
    if (((curve(node_index,2) <= geometry_torlance) && ...
            (curve(node_next_index,2) >= geometry_torlance)) || ...
            ((curve(node_index,2) >= geometry_torlance) && ...
            (curve(node_next_index,2) <= geometry_torlance)) )
        cross_flag=1;
        break;
    else
        cross_flag=0;
    end
end
end

function cross_flag=judgeCrossY(curve,geometry_torlance)
% function to judge if curve cross Y
%
cross_flag=0;
node_number=size(curve,1);
for node_index=1:node_number
    node_next_index=node_index+1;
    if node_next_index > node_number
        node_next_index=1;
    end
    if (((curve(node_index,1) <= geometry_torlance) && ...
            (curve(node_next_index,1) >= geometry_torlance)) || ...
            ((curve(node_index,1) >= geometry_torlance) && ...
            (curve(node_next_index,1) <= geometry_torlance)) )
        cross_flag=1;
        break;
    else
        cross_flag=0;
    end
end
end

function drawFlowField...
    (bias,jacobian)

a1=bias(1);
b1=jacobian(1);
c1=jacobian(3);
a2=bias(2);
b2=jacobian(2);
c2=jacobian(4);

x_list=(-1:0.1:1);
y_list=(-1:0.1:1);
hold on;
for x_index=1:length(x_list)
    x=x_list(x_index);
    for y_index=1:length(y_list)
        y=y_list(y_index);
        position=[x,y];
        vector=[a1+b1*x+c1*y,a2+b2*x+c2*y];
        quiver(position(1),position(2),vector(1),vector(2),'AutoScaleFactor',0.1);
    end
end
hold off

axis equal;
end

function drawFlowFieldProj...
    (bias,jacobian,stagnation_point,eig_vector)

a1=bias(1);
b1=jacobian(1);
c1=jacobian(3);
a2=bias(2);
b2=jacobian(2);
c2=jacobian(4);

x_list=(-1:0.1:1);
y_list=(-1:0.1:1);
hold on;
for x_index=1:length(x_list)
    x=x_list(x_index);
    for y_index=1:length(y_list)
        y=y_list(y_index);
        position=[x,y];
        vector=[a1+b1*x+c1*y,a2+b2*x+c2*y];
        position=(position-stagnation_point')/eig_vector';
        vector=vector/eig_vector';
        quiver(position(1),position(2),vector(1),vector(2),'AutoScaleFactor',0.1);
    end
end
hold off

axis equal;
end
