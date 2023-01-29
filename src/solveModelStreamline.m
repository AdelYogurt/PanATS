function solveModelStreamline()
% calculate element reference streaamline length
%
% reference:
% [1]KENWRIGHT D N. Automatic detection of open and closed
% separation and attachment lines; proceedings of the Proceedings
% Visualization '98 (Cat No98CB36276), F 18-23 Oct. 1998, 1998 [C].
%
global user_model

dimension=user_model.dimension;
point_list=user_model.point_list;
marker_list=user_model.marker_list;
MARKER_MONITERING=user_model.MARKER_MONITORING;

geometry_torlance=1e-12;

% calculate inflow vector
free_flow_vector=[1;0;0];

AOA=user_model.AOA/180*pi;
cos_AOA=cos(AOA);
sin_AOA=sin(AOA);
rotation_AOA=[
    cos_AOA 0 -sin_AOA;
    0 1 0;
    sin_AOA 0 cos_AOA];

AOS=user_model.SIDESLIP_ANGLE/180*pi;
cos_AOS=cos(AOS);
sin_AOS=sin(AOS);
rotation_SIDESLIP_ANGLE=[
    cos_AOS -sin_AOS 0;
    sin_AOS cos_AOS 0;
    0 0 1];

free_flow_vector=rotation_AOA*rotation_SIDESLIP_ANGLE*free_flow_vector;
user_model.free_flow_vector=free_flow_vector;

switch user_model.SYMMETRY
    case 'XOY'
        flip_operate=[1,1,-1];
        identify_dimension=3;
    case 'YOZ'
        flip_operate=[-1,1,1];
        identify_dimension=1;
    case 'ZOX'
        flip_operate=[1,-1,1];
        identify_dimension=2;
    otherwise
        error('solveModelStreamline: nuknown SYMMETRY type');
end

% calculate surface velocity
surface_flow_list=[zeros(size(point_list,1),4)]; % data, repeat times
for moniter_index=1:length(MARKER_MONITERING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITERING{moniter_index},marker_list);
    for element_index=1:marker_list(marker_index).element_number
        element=marker_element(element_index);
        normal_vector=element.normal_vector;
        point_index_list=element.point_index_list;
        dot_nor_vec=dot(free_flow_vector,normal_vector);
        Ve=(free_flow_vector'-normal_vector*dot_nor_vec);
        element.surface_flow=Ve;

        %         if judgeMatExistNum(point_index_list,92)
        %             disp('');
        %         end

        % add data to each point
        surface_flow_list(point_index_list,1:3)=...
            surface_flow_list(point_index_list,1:3)+Ve;
        surface_flow_list(point_index_list,4)=...
            surface_flow_list(point_index_list,4)+1;

        % identify which point of element on symmetry face
        % add symmetry velocity to point
        % (velocity identify_dimension set to zero)
        for node_index=1:length(point_index_list)
            point_index=point_index_list(node_index);
            if abs(point_list(point_index,identify_dimension)) <= ...
                    geometry_torlance
                % on symmetry face
                % add data to each point
                surface_flow_list(point_index,identify_dimension)=0;
            end
        end

        % process symmetry
        if element.element_nearby_number < length(element.point_index_list)
            % means element on symmetry face
            element_symmetry=HATSElement(element.element_type,...
                element.point_index_list);
            element_symmetry.center_point=element.center_point.*...
                flip_operate;
            element_symmetry.surface_flow=element.surface_flow.*...
                flip_operate;
            element_symmetry.normal_vector=element.normal_vector.*...
                flip_operate;

            % add to element_nearby_list
            element.element_nearby_list(element.element_nearby_number)=...
                element_symmetry;
        end
    end
end
% average velocity to each point
surface_flow_list(point_index_list,1:3)=...
    surface_flow_list(point_index_list,1:3)./surface_flow_list(point_index_list,4);

% record
user_model.surface_flow_list=surface_flow_list;

% identify stagnation element
% base on reference [1] to extract
for moniter_index=1:length(MARKER_MONITERING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITERING{moniter_index},marker_list);
    for element_index=1:marker_list(marker_index).element_number
        % base on around element to determine if stagnation element
        element=marker_element(element_index);
        if element.normal_vector*free_flow_vector > 0
            element.stagnation=int8(0);
            continue;
        end
        abs_Ve=norm(element.surface_flow);
        ez=element.normal_vector;
        if abs(abs_Ve) == 0
            element.stagnation=int8(0);
            continue;
        end

        % only consider upwind element have surface flow
        ex=element.surface_flow/abs_Ve;
        ey=cross(ez,ex);

        point_index_list=element.point_index_list;
        node_number=length(point_index_list);

        if element.normal_vector*free_flow_vector < -0.98
            %             close all hidden;
            %             surface_flow_node=surface_flow_list(point_index_list,1:3);
            %             drawLine(point_list(point_index_list,:))
            %             hold on;
            %             for node_index=1:node_number
            %                 point_index=point_index_list(node_index);
            %                 drawFlow(surface_flow_node(node_index,:),...
            %                     point_list(point_index,:))
            %             end
            %             hold off;
            %             view(3);
            %             xlabel('x');
            %             ylabel('y');
            %             zlabel('z');
            %             axis equal;
            %             disp('head');
        end

        % transform local cooridinate to 2D cooridinate
        % all vector just simply project to element panel
        Ve_list=zeros(node_number,3);
        d_node_list=zeros(node_number,3);
        for node_index=1:node_number
            point_index=point_index_list(node_index);
            Ve_list(node_index,:)=surface_flow_list(point_index,1:3)*[ex;ey;ez]';
            d_node_list(node_index,:)=(point_list(point_index,:)-...
                element.center_point)*[ex;ey;ez]';
        end
        d_node_list(:,end)=1;

        % interpolation fit flow
        coefficient=(d_node_list\Ve_list(:,[1,2]))';

        jacobian=coefficient(:,1:2); % [b1,c1;b2,c2]
        bias=coefficient(:,3); % [a1;a2]

        % jacobian matrix is negative, process stop
        det_jacobian=(jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1));
        if (((jacobian(1,1)+jacobian(2,2))^2-...
                4*det_jacobian) < 0)
            element.stagnation=int8(0);
            continue;
        end

        % evaluate eigenvalues of jacobian matrix
        % if one of eigenvalues is zero, process stop
        eig_value=eig(jacobian);
        if sum(abs(eig_value) <= geometry_torlance)
            element.stagnation=int8(0);
            continue;
        end

        if (element.normal_vector*free_flow_vector) < 0.98
%             drawElement...
%                 (d_node_list(:,[1,2]),Ve_list(:,[1,2]))
%             drawFlowField...
%                 (bias,jacobian);
%             disp('done');
        end

        stagnation_point=-jacobian\bias;
        % project element to canonical coordinates
        [eig_vector,eig_value]=eig(jacobian);
        eig_value=[eig_value(1,1),eig_value(2,2)];
        if (eig_value(1) < eig_value(2))
            eig_value=fliplr(eig_value);
            eig_vector=fliplr(eig_vector);
        end
        % normalize eig_vector
        d_node_list=(d_node_list(:,1:2)-stagnation_point')/eig_vector';

        if (element.normal_vector*free_flow_vector) < 0.98
%             drawElementProj...
%                 (d_node_list(:,[1,2]),Ve_list(:,[1,2]),stagnation_point,eig_vector)
%             drawFlowFieldProj...
%                 (bias,jacobian,stagnation_point,eig_vector);
%             disp('done');
        end

        % if (0,0) inside element
        surround_flag=judgeOriginSurround(d_node_list(:,[1,2]),geometry_torlance);
        if surround_flag
            element.stagnation=int8(1);
            continue;
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
            if judgeCrossY(d_node_list,node_number,geometry_torlance)
                element.stagnation=int8(1);
            else
                element.stagnation=int8(0);
            end
        end
    end
end

% % calculate streamline length
% inflow_direction_list=zeros(marker_element_number,3);
% streamline_length_list=zeros(marker_element_number,1);
%
% bottom_index_list=[];
% for element_index=1:marker_element_number
%     center_point=ADtree_marker_element.center_point_list(element_index,:);
%     out_index_list=HATS_element_list(element_index).out_index_list;
%     normal_vector=HATS_element_list(element_index).normal_vector;
%
%     if abs(normal_vector*free_flow_vector-1) < geometry_torlance
%         bottom_index_list=[bottom_index_list,element_index];
%         continue;
%     end
%
%     inflow_direction=inflow_direction_list(element_index,:);
%     streamline_length=streamline_length_list(element_index);
%
%     for out_node_index__=1:length(out_index_list)
%         out_node_index=out_index_list(out_node_index__);
%
%         out_center_point=ADtree_marker_element.center_point_list(out_node_index,:);
%
%         % base on normal_vector calculate inflow_direction and streamline_length
%         inflow_direction_out=(out_center_point-center_point);
%         streamline_length_out=sqrt(sum(inflow_direction_out.^2));
%         inflow_direction_out=inflow_direction_out/streamline_length_out;
%
%         % correct length and vector base on parent inflow_direction and streamline_length
%         if ~isempty(inflow_direction)
%             inflow_direction_out=...
%                 inflow_direction_out*streamline_length_out+...
%                 inflow_direction*streamline_length; % vector plus
%             streamline_length_out=sqrt(sum(inflow_direction_out.^2));
%             inflow_direction_out=inflow_direction_out/streamline_length_out;
%         end
%
%         % compare exist steamline, if shorter than repace
%         if (streamline_length_list(out_node_index) ~= 0)
%             if streamline_length_list(out_node_index) > streamline_length_out
%                 inflow_direction_list(out_node_index,:)=inflow_direction_out;
%                 streamline_length_list(out_node_index)=streamline_length_out;
%             end
%         else
%             inflow_direction_list(out_node_index,:)=inflow_direction_out;
%             streamline_length_list(out_node_index)=streamline_length_out;
%         end
%     end
% end
% inflow_direction_list(bottom_index_list,:)=repmat(free_flow_vector',length(bottom_index_list),1);
% streamline_length_list(bottom_index_list)=max(streamline_length_list);
%
% streamline_output.inflow_direction_list=inflow_direction_list;
% streamline_output.streamline_length_list=streamline_length_list;
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

function drawElement(curve,vector_list)
node_number=size(curve,1);
order=[1:node_number,1];
line(curve(order,1),curve(order,2));
hold on;
for node_index=1:node_number
    position=curve(node_index,:);
    vector=vector_list(node_index,:);
    quiver(position(1),position(2),vector(1),vector(2),'AutoScaleFactor',0.001);
end
hold off;
end

function drawElementProj(curve,vector_list,stagnation_point,eig_vector)
node_number=size(curve,1);
order=[1:node_number,1];
line(curve(order,1),curve(order,2));
hold on;
for node_index=1:node_number
    position=curve(node_index,:);
    vector=vector_list(node_index,:);
    vector=vector/eig_vector';
    quiver(position(1),position(2),vector(1),vector(2),'AutoScaleFactor',0.001);
end
hold off;
end

function drawFlowField...
    (bias,jacobian)

a1=bias(1);
b1=jacobian(1);
c1=jacobian(3);
a2=bias(2);
b2=jacobian(2);
c2=jacobian(4);

x_list=(-1:0.1:1)*1e-2;
y_list=(-1:0.1:1)*1e-2;
hold on;
for x_index=1:length(x_list)
    x=x_list(x_index);
    for y_index=1:length(y_list)
        y=y_list(y_index);
        position=[x,y];
        vector=[a1+b1*x+c1*y,a2+b2*x+c2*y];
        quiver(position(1),position(2),vector(1),vector(2),'AutoScaleFactor',0.0001);
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

x_list=(-1:0.1:1)*1e-2;
y_list=(-1:0.1:1)*1e-2;
hold on;
for x_index=1:length(x_list)
    x=x_list(x_index);
    for y_index=1:length(y_list)
        y=y_list(y_index);
        position=[x,y];
        vector=[a1+b1*x+c1*y,a2+b2*x+c2*y];
        position=(position-stagnation_point')/eig_vector';
        vector=vector/eig_vector';
        quiver(position(1),position(2),vector(1),vector(2),'AutoScaleFactor',0.0001);
    end
end
hold off

axis equal;
end
