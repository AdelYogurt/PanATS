clc;
clear;
close all hidden;

geometry_torlance=1e-12;

order=[1,2,3,1];

element=HATSElement(5,[1,2,3]);

d_node_list=[
    -1,-1,1;
    1,-1,1;
    0,1,1];
Ve_list=[
    1,-1;
    -1,-1;
    0,2];
node_number=3;

% interpolation fit flow
coefficient=(d_node_list\Ve_list(:,[1,2]))';

jacobian=coefficient(:,1:2); % [b1,c1;b2,c2]
bias=coefficient(:,3); % [a1;a2]

% jacobian matrix is negative, process stop
det_jacobian=(jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1));
if (((jacobian(1,1)+jacobian(2,2))^2-...
        4*det_jacobian) < 0)
    element.stagnation=int8(0);
end

% evaluate eigenvalues of jacobian matrix
% if one of eigenvalues is zero, process stop
eig_value=eig(jacobian);
if sum(abs(eig_value) <= geometry_torlance)
    element.stagnation=int8(0);
end

stagnation_point=-jacobian\bias;
% project element to canonical coordinates
[eig_vector,eig_value]=eig(jacobian);
eig_value=[eig_value(1,1),eig_value(2,2)];
if (eig_value(1) < eig_value(2))
    eig_value=fliplr(eig_value);
    eig_vector=fliplr(eig_vector);
end
d_node_list=(d_node_list(:,1:2)-stagnation_point')/eig_vector';

line(d_node_list(order,1),d_node_list(order,2));
drawFlowFieldProj...
    (bias,jacobian,stagnation_point,eig_vector);

% if (0,0) inside element
surround_flag=judgeOriginSurround(d_node_list(:,[1,2]),geometry_torlance);
if surround_flag
    element.stagnation=int8(1);
end

% judge phase portrait
% eigenvalue less than zero is concentrate
% X corresponds to the eigenvectors 1
% Y corresponds to the eigenvectors 2
if ((eig_value(1) > geometry_torlance) && ...
        (eig_value(2) > geometry_torlance))
    % repelling node, check small eigenvalue corresponded axis
    % if cross Y
    if judgeCrossY(d_node_list,node_number,geometry_torlance)
        element.stagnation=int8(1);
    else
        element.stagnation=int8(0);
    end
elseif ((eig_value(1) < geometry_torlance) && ...
        (eig_value(2) < geometry_torlance))
    % attracting node, check
    % concentrate line is not stagnation point
    element.stagnation=int8(0);
else
    % saddle, judge axis X and axis Y which is separation line
    % means which eigenvalue is large than zero
    % repelling node, check small eigenvalue corresponded axis
    % if cross X
    if judgeCrossY(d_node_list,node_number,geometry_torlance)
        element.stagnation=int8(1);
    else
        element.stagnation=int8(0);
    end
end

function surround_flag=judgeOriginSurround(curve,geometry_torlance)
% function to judge if origin point(0,0) is surrounded by line
%
surround_flag=0;
cross_time=0;
for edge_index=1:size(curve,1)
    edge_index_next=edge_index+1;
    if edge_index_next > size(curve,1)
        edge_index_next=1;
    end

    if ((-curve(edge_index_next,2)) * ...
            ((curve(edge_index,2)))) < -eps
        % if cross point not in line y range
        continue;
    end

    A=curve(edge_index,2)-curve(edge_index_next,2);
    C=curve(edge_index_next,2)*curve(edge_index,1)-...
        curve(edge_index_next,1)*curve(edge_index,2);
    cross_point=-C/A;
    if abs(A) > geometry_torlance
        if cross_point > 0
            if ((cross_point - curve(edge_index_next,1)) * ...
                    ((curve(edge_index,1)) - cross_point)) >= -geometry_torlance
                % cross point must in line range
                cross_time=cross_time+1;
            end
        end
    else
        if ((-curve(edge_index_next,1)) * ...
                ((curve(edge_index,1)))) >= geometry_torlance
            % cross point must in line range
            % line overlap X axis
            cross_time=cross_time+1;
        end
    end
end
if (cross_time == 1)
    surround_flag=1;
end
end

function cross_flag=judgeCrossX(curve,node_number,geometry_torlance)
% function to judge if curve cross X
%
cross_flag=0;
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

function cross_flag=judgeCrossY(curve,node_number,geometry_torlance)
% function to judge if curve cross Y
%
cross_flag=0;
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

