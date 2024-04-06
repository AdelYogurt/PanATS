function [cross_flag,stagnation_flag]=judgeAttachment(normal_vector,node_list,node_flow_list)
% judge if this element have attachment line across
% 
% reference: [1] KENWRIGHT D N. Automatic detection of open and closed
% separation and attachment lines; proceedings of the Proceedings
% Visualization '98 (Cat No98CB36276), F 18-23 Oct. 1998, 1998 [C].
%
precision_torlance=eps;
topology_torlance=1e-2;

cross_flag=false(1);
stagnation_flag=false(1);

% project to 2D
e3=normal_vector;
e1=node_list(1,1:3)-node_list(2,1:3);
e1=e1-(e1*e3')*e3; % orthogonalization
e1=e1/norm(e1);
e2=cross(e3,e1);

% transform local cooridinate to 2D cooridinate
% project flow to plane, and scale norm of flow to initial flow
Len_flow_init=vecnorm(node_flow_list,2,2);
node_flow_list=node_flow_list*[e1;e2;e3]';
Len_flow_proj=vecnorm(node_flow_list,2,2);Len_flow_proj(Len_flow_proj == 0)=1;
node_flow_list=node_flow_list./Len_flow_proj.*Len_flow_init;
node_list=node_list*[e1;e2;e3]';
node_list(:,end)=1;

% interpolation fit flow
coefficient=(node_list\node_flow_list(:,[1,2]))';

jacobian=coefficient(:,1:2); % [b1,c1;b2,c2]
% jacobian(abs(jacobian) < precision_torlance)=0; % solve precision, avoid negtive
bias=coefficient(:,3); % [a1;a2]

% % velocity do not too close
% low_bou=min(node_list(:,1:2),[],1);
% up_bou=max(node_list(:,1:2),[],1);
% center_bou=(low_bou+up_bou)/2;
% aver_bou=sum(up_bou-low_bou)/2;
% point_1=center_bou+[aver_bou,aver_bou];
% point_2=center_bou+[-aver_bou,aver_bou];
% Ve_1=point_1*jacobian'+bias';
% Ve_2=point_2*jacobian'+bias';
% cos_angle=Ve_1*Ve_2'/norm(Ve_1)/norm(Ve_2);
% if cos_angle > 0.8
%     return;
% end

% jacobian matrix is negative, phase is spiral
det_jacobian=(jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1));
if (((jacobian(1,1)+jacobian(2,2))^2-4*det_jacobian) < -precision_torlance)
    % check if spiral is out flow
    eig_value=eig(jacobian);
    if real(eig_value(1)) > 0
        % out flow
        stagnation_point=-jacobian\bias;
        node_proj_list=(node_list(:,1:2)-stagnation_point');
        offset = norm(max(node_proj_list)-min(node_proj_list))*topology_torlance;
        node_proj_list=geoCurveOffset(node_proj_list,offset);

        % if (0,0) inside element
        if judgeOriginSurround(node_proj_list)
            cross_flag=true(1);
            stagnation_flag=true(1);
        end
    end
    return;
end

% project element to canonical coordinates
[eig_vector,eig_value]=eig(jacobian);
eig_value=[eig_value(1,1),eig_value(2,2)];
if (eig_value(1) < eig_value(2))
    eig_value=fliplr(eig_value);
    eig_vector=fliplr(eig_vector);
end
% from now, eig_value(1) > eig_value(2)
% X corresponds to the eigenvectors 1
% Y corresponds to the eigenvectors 2

% evaluate eigenvalues of jacobian matrix
% if one of eigenvalues is zero, process stop
if all(abs(eig_value) < topology_torlance)
    stagnation_point=[0;0];
    return;
elseif abs(eig_value(2)) < precision_torlance
    % check if cross separation line
    if abs(jacobian(1)) < precision_torlance
        x_c=-bias(1);
    else
        x_c=-bias(1)/jacobian(1);
    end
    if abs(jacobian(4)) < precision_torlance
        y_c=-bias(2);
    else
        y_c=-bias(2)/jacobian(4);
    end
    stagnation_point=[x_c;y_c];
elseif abs(eig_value(1)) < precision_torlance
    % converge line
    if abs(jacobian(1)) < precision_torlance
        x_c=-bias(1);
    else
        x_c=-bias(1)/jacobian(1);
    end
    if abs(jacobian(4)) < precision_torlance
        y_c=-bias(2);
    else
        y_c=-bias(2)/jacobian(4);
    end
    stagnation_point=[x_c;y_c];
    return;
else
    if rcond(jacobian) < eps
        stagnation_point=lsqminnorm(jacobian,-bias);
    else
        stagnation_point=-jacobian\bias;
    end
end

% normalize point_2D
node_proj_list=(node_list(:,1:2)-stagnation_point')/eig_vector';

% velocity do not too close
max_bou=max(node_proj_list,[],1);
min_bou=min(node_proj_list,[],1);
center_bou=(min_bou+max_bou)/2;
aver_bou=sum(max_bou-min_bou)/2;
point_1=center_bou+[aver_bou,aver_bou];
point_2=center_bou+[-aver_bou,aver_bou];
Ve_1=eig_value.*point_1;
Ve_2=eig_value.*point_2;
cos_angle=Ve_1*Ve_2'/norm(Ve_1)/norm(Ve_2);
if cos_angle > 0.8
    return;
end

% % normalize, becase amplitude now is useless
% node_proj_list(:,1)=node_proj_list(:,1)/max(abs(node_proj_list(:,1)));
% node_proj_list(:,2)=node_proj_list(:,2)/max(abs(node_proj_list(:,2)));

bou_node=node_proj_list(1:end-1,:);
if det(eig_vector) < 0 % if det(eig_vector) less than 0, rotation
    bou_node=flipud(bou_node);
end
offset=norm(max(bou_node)-min(bou_node))*topology_torlance;
bou_node=geoCurveOffset(bou_node,offset);

% judge phase portrait
% eigenvalue less than zero is concentrate
% X corresponds to the eigenvectors 1
% Y corresponds to the eigenvectors 2
if ((eig_value(1) > 0) && (eig_value(2) >= 0))
    if abs(eig_value(1)-eig_value(2)) < topology_torlance
        % proper node
        % if (0,0) inside element
        if judgeOriginSurround(bou_node)
            cross_flag=true(1);
            stagnation_flag=true(1);
        end
    else
        % repelling node, check small eigenvalue corresponded axis
        % if cross Y
        if judgeCrossY(bou_node(1:end,[1,2]))
            cross_flag=true(1);
        end
    end
elseif ((eig_value(1) < -0) && (eig_value(2) < -0))
    % attracting node, check
    % concentrate line is not stagnation point
    
elseif ((eig_value(1) > 0) && (eig_value(2) < -0))
    % saddle, judge axis X and axis Y which is separation line
    % means which eigenvalue is large than zero
    % repelling node, check small eigenvalue corresponded axis
    if judgeCrossY(bou_node(1:end,[1,2]))
        cross_flag=true(1);
    end
end

end

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

    if norm(dr) == 0
        return;
    end

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

function surround_flag=judgeOriginSurround(curve)
% function to judge if origin point(0,0) is surrounded by line
% calculate curve cross positive axis x times, if is odd number, surround
%
surround_flag=0;

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

function cross_flag=judgeCrossX(curve)
% function to judge if curve cross X
%
cross_flag=false(1);
node_number=size(curve,1);
for node_index=1:node_number
    node_next_index=node_index+1;
    if node_next_index > node_number
        node_next_index=1;
    end
    if (((curve(node_index,2) <= 0) && ...
            (curve(node_next_index,2) >= 0)) || ...
            ((curve(node_index,2) >= 0) && ...
            (curve(node_next_index,2) <= 0)) )
        cross_flag=true(1);
        break;
    else
        cross_flag=false(1);
    end
end
end

function cross_flag=judgeCrossY(curve)
% function to judge if curve cross Y
%
cross_flag=false(1);
node_number=size(curve,1);
for node_index=1:node_number
    node_next_index=node_index+1;
    if node_next_index > node_number
        node_next_index=1;
    end
    if (((curve(node_index,1) <= 0) && ...
            (curve(node_next_index,1) >= 0)) || ...
            ((curve(node_index,1) >= 0) && ...
            (curve(node_next_index,1) <= 0)) )
        cross_flag=true(1);
        break;
    else
        cross_flag=false(1);
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
