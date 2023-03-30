function solveModelStreamline()
% identify specifical element(stagnation elemenet, separate element)
% specifical element is the start of streamline
% calculate element reference streaamline length
%
% point_list is coordinate of all node
% element_list contain []
% marker_list contain make{marker_name,marker_element_number,marker_element}
% marker_element contain HATSElement
%
% reference: [1] KENWRIGHT D N. Automatic detection of open and closed
% separation and attachment lines; proceedings of the Proceedings
% Visualization '98 (Cat No98CB36276), F 18-23 Oct. 1998, 1998 [C].
%
% copyright Adel 2023.03
%
global user_model

dimension=user_model.dimension;

point_list=user_model.point_list;
edge_list=user_model.edge_list;
marker_list=user_model.marker_list;

MARKER_MONITERING=user_model.MARKER_MONITORING;

geometry_torlance=1e-9;

% calculate free flow vector
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

% clear last time record
for monitor_index=1:length(MARKER_MONITERING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITERING(monitor_index),marker_list);
    for element_index=1:marker_list(marker_index).element_number
        marker_element(element_index).streamline_ref=[];
    end
end

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

% calculate surface velocity base on velocity project
surface_flow_list=[zeros(size(point_list,1),4)]; % surface_flow, repeat times of vertex
for moniter_index=1:length(MARKER_MONITERING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITERING{moniter_index},marker_list);
    for element_index=1:marker_list(marker_index).element_number
        element=marker_element(element_index);
        normal_vector=element.normal_vector;
        point_index_list=element.point_index_list;
        element.attachment=[];

        % calculate surface flow of element(norm value less than 1)
        dot_nor_vec=dot(free_flow_vector,normal_vector);
        if abs(dot_nor_vec-1) < geometry_torlance
            Ve=[0,0,0];
        else
            Ve=(free_flow_vector'-normal_vector*dot_nor_vec);
        end

        element.surface_flow=Ve;

        % add data to each point
        surface_flow_list(point_index_list,1:3)=...
            surface_flow_list(point_index_list,1:3)+Ve;
        surface_flow_list(point_index_list,4)=...
            surface_flow_list(point_index_list,4)+1;

        % identify which point of element on symmetry face
        % add symmetry velocity to point
        % (velocity identify_dimension set to zero)
        for node_index=1:length(point_index_list)
            vertex_index=point_index_list(node_index);
            if abs(point_list(vertex_index,identify_dimension)) <= ...
                    geometry_torlance
                % on symmetry face
                % add data to each point
                surface_flow_list(vertex_index,identify_dimension)=0;
            end
        end

%         % process symmetry
%         if element.element_nearby_number < length(element.point_index_list)
%             % means element on symmetry face
%             element_symmetry=HATSElement(element.element_type,...
%                 element.point_index_list);
%             element_symmetry.center_point=element.center_point.*...
%                 flip_operate;
%             element_symmetry.surface_flow=element.surface_flow.*...
%                 flip_operate;
%             element_symmetry.normal_vector=element.normal_vector.*...
%                 flip_operate;
% 
%             % add to element_nearby_list
%             element.element_nearby_list(element.element_nearby_number+1)=...
%                 element_symmetry;
%         end
    end
end

boolean = surface_flow_list(:,4) ~= 0;
% average velocity to each point
surface_flow_list(boolean,1:3)=...
    surface_flow_list(boolean,1:3)./surface_flow_list(boolean,4);

% record
user_model.surface_flow_list=surface_flow_list;

% identify attachment element
attachment_list=HATSElement.empty();
for moniter_index=1:length(MARKER_MONITERING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITERING{moniter_index},marker_list);
    for element_index=1:marker_list(marker_index).element_number
        element=marker_element(element_index);
        normal_vector=element.normal_vector;
        point_arou_list=point_list(element.point_index_list,:);
        Ve_list=surface_flow_list(element.point_index_list,1:3);

        if dot(normal_vector,free_flow_vector) < 0 %% only upwind element
            [cross_flag,stagnation_flag]=judgeAttachment(normal_vector,point_arou_list,Ve_list);
            if cross_flag
                element.attachment=1;
                attachment_list=[attachment_list;element];
            else
                element.attachment=0;
            end

%             if stagnation_flag
%                 stagnation_list=[stagnation_list;moniter_index,element_index];
%             end
        else
            element.attachment=0;
        end
    end
end

% initialize data sort array
% streamline_len_list, streamline_list
output_streamline.streamline_len_list=cell(length(marker_list),1);
output_streamline.streamline_list={};
output_streamline.attachment_list=attachment_list;

max_streamline_len=0; % max length of streamline

% select stagnation element center point as start point to generate streamline
for attachment_index=1:length(attachment_list)
    element=attachment_list(attachment_index);

    streamline_index=length(output_streamline.streamline_list)+1;
    streamline=calStreamline...
        (free_flow_vector,point_list,edge_list,element,element.center_point,...
        user_model.element_number,streamline_index,geometry_torlance);
%     line(streamline(:,1),streamline(:,2),streamline(:,3),'Color','r');

    output_streamline.streamline_list{streamline_index,1}=streamline;
    if streamline(end,4) > max_streamline_len
        max_streamline_len=streamline(end,4);
    end
end

% check which element do not cross by streamline
% if exist, add streamline that cross this element
for monitor_index=1:length(MARKER_MONITERING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITERING(monitor_index),marker_list);

    for element_index=1:marker_list(marker_index).element_number
        element=marker_element(element_index);

        if isempty(element.streamline_ref) && ...
                sum(abs(element.surface_flow)) > geometry_torlance
            % if donot have streamline cross, add new streamline
            streamline_index=length(output_streamline.streamline_list)+1;
            streamline=calStreamline...
                (free_flow_vector,point_list,edge_list,element,element.center_point,...
                user_model.element_number,streamline_index,geometry_torlance);
%             line(streamline(:,1),streamline(:,2),streamline(:,3),'Color','r');

            output_streamline.streamline_list{streamline_index,1}=streamline;
            if streamline(end,4) > max_streamline_len
                max_streamline_len=streamline(end,4);
            end

%             if streamline_index == 72
%                 disp('?');
%             end
        end
    end
end

for monitor_index=1:length(MARKER_MONITERING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITERING(monitor_index),marker_list);

    streamline_len_list=zeros(marker_list(marker_index).element_number,1);

    for element_index=1:marker_list(marker_index).element_number
        element=marker_element(element_index);
        streamline_ref=element.streamline_ref;

        if isempty(streamline_ref)
            if sum(abs(element.surface_flow)) < geometry_torlance
                % for backup element and have no surface_flow, set to max_streamline_len
                streamline_len_list(element_index)=max_streamline_len;
            else
                disp('?');
            end
        else
            streamline_ref_number=size(streamline_ref,1);
            if streamline_ref_number == 1
                streamline=output_streamline.streamline_list{streamline_ref(1)};
                base_len=streamline(streamline_ref(2),4);
                base_point=streamline(streamline_ref(2),1:3);

                streamline_len=base_len+norm(element.center_point-base_point);
            else
                % select min streamline length as streamline_ref
                streamline_len=max_streamline_len;
                for streamline_ref_index=1:streamline_ref_number
                    streamline=output_streamline.streamline_list{streamline_ref(streamline_ref_index,1)};
                    base_len=streamline(streamline_ref(streamline_ref_index,2),4);
                    base_point=streamline(streamline_ref(streamline_ref_index,2),1:3);

                    streamline_len_temp=base_len+norm(element.center_point-base_point);
                    if streamline_len_temp < streamline_len
                        streamline_len=streamline_len_temp;
                    end
                end
            end

            streamline_len_list(element_index)=streamline_len;
        end
    end

    output_streamline.streamline_len_list{marker_index}=streamline_len_list;
end

user_model.output_streamline=output_streamline;

if user_model.INFORMATION
    fprintf('solveModelStreamline: streamline solve done!\n');
    fprintf('max streamline length: %14f\n',max_streamline_len);
end

end

function streamline=calStreamline...
    (free_flow_vector,point_list,edge_list,element,point_start,...
    max_size,streamline_index,geometry_torlance)
% calculate streamline from point_start and element
%
streamline=zeros(max_size,4); % cross point, current length
edge_cross_list=zeros(max_size,2,'int32'); % cross edge vertex, use for reverse search
index=1;
streamline(index,:)=[point_start,0];

% upstream search begin

% if element is no attachment, begin upstream search
if ~element.attachment
    surface_flow=element.surface_flow;
    point_index_list=element.point_index_list;
    point_number=length(point_index_list);

    % edge should place on upstream of point_start
    % downstream edge will become start edge
    for point_index=1:point_number
        vertex_index=point_index_list(point_index);
        if point_index == point_number
            vertex_ref_index=point_index_list(1);
        else
            vertex_ref_index=point_index_list(point_index+1);
        end
        point=point_list(vertex_index,1:3);
        point_ref=point_list(vertex_ref_index,1:3);

        dr=point-point_start;
        dr_ref=point_ref-point_start;
        if ((dr/norm(dr)+dr_ref/norm(dr_ref))*surface_flow') >= 0
            if cross(surface_flow,dr)*cross(surface_flow,dr_ref)' <= 0 % edge should cross surface_flow
                vertex_forbid_index=vertex_index;
                break;
            end
        end
    end

    upstream_flag=true(1);
end

% upstream cross point search
while ~isempty(element) && index < max_size && ...
        ~element.attachment && upstream_flag
    % load element information
    point_index_list=element.point_index_list;
    surface_flow=element.surface_flow;
    point_number=length(point_index_list);

    % calculate next cross point
    % calculate streamline cross point on edge
    % vertex_index, vertex_ref_index is edge than this surface_flow cross
    point_end=[];
    for point_index=1:point_number
        vertex_index=point_index_list(point_index);

        % start edge donot need to check
        if vertex_index == vertex_forbid_index
            continue;
        end

        if point_index == point_number
            vertex_ref_index=point_index_list(1);
        else
            vertex_ref_index=point_index_list(point_index+1);
        end
        point=point_list(vertex_index,1:3);
        point_ref=point_list(vertex_ref_index,1:3);

        dr=point-point_start;
        dr_ref=point_ref-point_start;

        % check if edge was upstream 
        if ((dr/norm(dr)+dr_ref/norm(dr_ref))*surface_flow') <= 0
            % check if edge was crossed by surface_flow
            if cross(surface_flow,dr)*cross(surface_flow,dr_ref)' <= 0
                point_end=calCrossPoint(point,point_ref,point_start,surface_flow);
                break;
            end
        end
    end

    if isempty(point_end)
        % means cannot upstream from face, search by edge
        dedge_ref_list=(point_list(point_index_list,1:3)-point_start);
        dedge_ref_list=dedge_ref_list./sqrt(sum(dedge_ref_list.^2,2));
        proj_vec=dedge_ref_list*free_flow_vector;
        [~,upstream_index]=min(proj_vec);
        vertex_index=point_index_list(upstream_index);
        if upstream_index == length(point_index_list)
            vertex_ref_index=point_index_list(1);
        else
            vertex_ref_index=point_index_list(upstream_index+1);
        end
        point_end=point_list(vertex_index,1:3);
        %     line(point_end(:,1),point_end(:,2),point_end(:,3),'Marker','o')

        % record parameter
        index=index+1;
        streamline(index,:)=[point_end,streamline(index-1,4)-norm(point_start-point_end)];
        edge_cross_list(index-1,:)=[vertex_index,vertex_ref_index];

        % search next element
        edge=edge_list(vertex_index);
        
        % search if exist element that surface flow is between two edge
        element_search_flag=false(1);
        for edge_index=1:edge.edge_number
            element_ref=edge.element_ref_list(edge_index);
            vertex_ref_index=edge.vertex_ref_list(edge_index);
            [~,point_index]=judgeMatExistNum...
                (element_ref.point_index_list,vertex_index);
            if point_index == 1
                vertex_prev_index=element_ref.point_index_list(end);
            else
                vertex_prev_index=element_ref.point_index_list(point_index-1);
            end

            dedge_ref=point_list(vertex_ref_index,1:3)-point_end;
            dedge_prev=point_list(vertex_prev_index,1:3)-point_end;

            % check if edge was upstream 
            if ((dedge_ref/norm(dedge_ref)+dedge_prev/norm(dedge_prev))*element_ref.surface_flow') <= 0
                % check if edge was crossed by surface_flow
                if cross(element_ref.surface_flow,dedge_ref)*cross(element_ref.surface_flow,dedge_prev)' <= 0 % edge should cross surface_flow
                    element_search_flag=true(1);
                    break;
                end
            end
        end

        % if cannot find element that surface flow is between two edge
        % upstream search from edge
        if ~element_search_flag
            dedge_ref_list=point_list(edge.vertex_ref_list(1:edge.edge_number),1:3)-point_end;
            dedge_ref_list=dedge_ref_list./sqrt(sum(dedge_ref_list.^2,2));
            proj_edge=dedge_ref_list*free_flow_vector;
            [~,point_index]=min(proj_edge);
            element_ref=edge.element_ref_list(point_index);
        end

        % next search element
        element=element_ref;
        vertex_forbid_index=0;
        point_start=point_end;
    else
        %     line(point_end(:,1),point_end(:,2),point_end(:,3),'Marker','o')

        % record parameter
        index=index+1;
        streamline(index,:)=[point_end,streamline(index-1,4)-norm(point_start-point_end)];
        edge_cross_list(index-1,:)=[vertex_index,vertex_ref_index];

        % search element_ref of edge opposite
        % next search element
        edge=edge_list(vertex_ref_index);
        element_ref=edge.element_ref_list(edge.getRefIndex(vertex_index));

        if ~isempty(element_ref)
            element=element_ref;
        end

        vertex_forbid_index=vertex_ref_index;
        point_start=point_end;
    end
end

% reverse and connect point
if index > 1
    if isempty(element)
        % call back last time element
        edge=edge_list(edge_cross_list(index-1,1));
        element=edge.element_ref_list(edge.getRefIndex(edge_cross_list(index-1,2)));
        streamline(index,:)=zeros(1,4);
        edge_cross_list(index-1,:)=zeros(1,2,'int32');
        index=index-1;
    else
        if element.attachment
            % load element information
            point_index_list=element.point_index_list;
            surface_flow=element.surface_flow;
            point_number=length(point_index_list);

            % calculate next cross point
            point_end=[];
            for point_index=1:point_number
                vertex_index=point_index_list(point_index);

                % start edge donot need to check
                if vertex_index == vertex_forbid_index
                    continue;
                end

                if point_index == point_number
                    vertex_ref_index=point_index_list(1);
                else
                    vertex_ref_index=point_index_list(point_index+1);
                end
                point=point_list(vertex_index,1:3);
                point_ref=point_list(vertex_ref_index,1:3);

                dr=point-point_start;
                dr_ref=point_ref-point_start;
                if cross(surface_flow,dr)*cross(surface_flow,dr_ref)' <= 0 % edge should cross surface_flow
                    point_end=calCrossPoint(point,point_ref,point_start,surface_flow);
                    break;
                end
            end

            % record parameter
            index=index+1;
            streamline(index,:)=[point_end,streamline(index-1,4)-norm(point_start-point_end)];
            edge_cross_list(index-1,:)=[vertex_index,vertex_ref_index];
        end
    end

    max_length=-streamline(index,4);
    streamline(1:index,:)=streamline(index:-1:1,:);

    streamline(1:index,4)=streamline(1:index,4)+max_length;

    index=index-1;
    point_start=streamline(index,1:3);

    % add element.streamline_ref downstream
    for down_index=1:index
        element.streamline_ref=[element.streamline_ref;...
            streamline_index,down_index];
        last_index=index-down_index+1;
        edge=edge_list(edge_cross_list(last_index,1));
        element=edge.element_ref_list(edge.getRefIndex(edge_cross_list(last_index,2)));
    end

    vertex_forbid_index=edge_cross_list(last_index,1);
else
    % downstream search begin
    surface_flow=element.surface_flow;
    point_index_list=element.point_index_list;
    point_number=length(point_index_list);

    % edge should place on upstream of point_start
    % downstream edge will become start edge
    for point_index=1:point_number
        vertex_index=point_index_list(point_index);
        if point_index == point_number
            vertex_ref_index=point_index_list(1);
        else
            vertex_ref_index=point_index_list(point_index+1);
        end
        point=point_list(vertex_index,1:3);
        point_ref=point_list(vertex_ref_index,1:3);

        dr=point-point_start;
        dr_ref=point_ref-point_start;
        if ((dr/norm(dr)+dr_ref/norm(dr_ref))*surface_flow') <= 0
            if cross(surface_flow,dr)*cross(surface_flow,dr_ref)' <= 0
                vertex_forbid_index=vertex_index;
                break;
            end
        end
    end

end

% downstream search begin

downstream_flag=true(1);
while index < max_size && downstream_flag
    % load element information
    point_index_list=element.point_index_list;
    surface_flow=element.surface_flow;
    point_number=length(point_index_list);

    % calculate next cross point
    point_end=[];
    for point_index=1:point_number
        vertex_index=point_index_list(point_index);

        % start edge donot need to check
        if vertex_index == vertex_forbid_index
            continue;
        end

        if point_index == point_number
            vertex_ref_index=point_index_list(1);
        else
            vertex_ref_index=point_index_list(point_index+1);
        end
        point=point_list(vertex_index,1:3);
        point_ref=point_list(vertex_ref_index,1:3);

        dr=point-point_start;
        dr_ref=point_ref-point_start;

        % check if edge was downstream 
        if ((dr/norm(dr)+dr_ref/norm(dr_ref))*surface_flow') >= 0
            % check if edge was crossed by surface_flow
            if cross(surface_flow,dr)*cross(surface_flow,dr_ref)' <= 0
                point_end=calCrossPoint(point,point_ref,point_start,surface_flow);
                break;
            end
        end
    end

    if isempty(point_end)
        % means cannot downstream from face, search by edge
        dedge_ref_list=(point_list(point_index_list,1:3)-point_start);
        dedge_ref_list=dedge_ref_list./sqrt(sum(dedge_ref_list.^2,2));
        proj_vec=dedge_ref_list*free_flow_vector;
        [~,upstream_index]=max(proj_vec);
        vertex_index=point_index_list(upstream_index);
        if upstream_index == length(point_index_list)
            vertex_ref_index=point_index_list(1);
        else
            vertex_ref_index=point_index_list(upstream_index+1);
        end
        point_end=point_list(vertex_index,1:3);
        %     line(point_end(:,1),point_end(:,2),point_end(:,3),'Marker','o')

        % if search point is reverse flow, end search
        if (point_end-point_start)*free_flow_vector <= geometry_torlance
            break;
        end

        % record parameter
        index=index+1;
        streamline(index,:)=[point_end,streamline(index-1,4)+norm(point_start-point_end)];
        edge_cross_list(index-1,:)=[vertex_index,vertex_ref_index];

        % search next element
        edge=edge_list(vertex_index);
        
        % search if exist element that surface flow is between two edge
        element_search_flag=false(1);
        for edge_index=1:edge.edge_number
            element_ref=edge.element_ref_list(edge_index);
            vertex_ref_index=edge.vertex_ref_list(edge_index);
            [~,point_index]=judgeMatExistNum...
                (element_ref.point_index_list,vertex_index);
            if point_index == 1
                vertex_prev_index=element_ref.point_index_list(end);
            else
                vertex_prev_index=element_ref.point_index_list(point_index-1);
            end

            dedge_ref=point_list(vertex_ref_index,1:3)-point_end;
            dedge_prev=point_list(vertex_prev_index,1:3)-point_end;

            % check if edge was downstream 
            if ((dedge_ref/norm(dedge_ref)+dedge_prev/norm(dedge_prev))*element_ref.surface_flow') >= 0
                % check if edge was crossed by surface_flow
                if cross(element_ref.surface_flow,dedge_ref)*cross(element_ref.surface_flow,dedge_prev)' <= 0
                    element_search_flag=true(1);
                    break;
                end
            end
        end

        % if cannot find element that surface flow is between two edge
        % downstream search from edge
        if ~element_search_flag
            dedge_ref_list=point_list(edge.vertex_ref_list(1:edge.edge_number),1:3)-point_end;
            dedge_ref_list=dedge_ref_list./sqrt(sum(dedge_ref_list.^2,2));
            proj_edge=dedge_ref_list*free_flow_vector;
            [~,point_index]=max(proj_edge);
            element_ref=edge.element_ref_list(point_index);
        end

    else
        %     line(point_end(:,1),point_end(:,2),point_end(:,3),'Marker','o')

        % record parameter
        element.streamline_ref=[element.streamline_ref;...
            streamline_index,index];
        index=index+1;
        streamline(index,:)=[point_end,streamline(index-1,4)+norm(point_start-point_end)];

        % search element_ref of edge opposite
        edge=edge_list(vertex_ref_index);
        element_ref=edge.element_ref_list(edge.getRefIndex(vertex_index));
    end

    if isempty(element_ref)
        downstream_flag=false(1);
    end

    % next iteration
    element=element_ref;
    vertex_forbid_index=vertex_ref_index;
    point_start=point_end;

end

streamline=streamline(1:index,:);
end

function point_end=calCrossPoint(point,point_ref,point_start,surface_flow)
% calculate cross point coordinate as surface_flow downstream point_end
%
surface_flow=surface_flow/norm(surface_flow);
dr_edge=point_ref-point;
dr_edge_length=norm(dr_edge);
dr_edge_norm=dr_edge/dr_edge_length;
point_proj_length=norm(cross(surface_flow,point-point_start));
point_end=point+point_proj_length/norm(cross(surface_flow,dr_edge_norm))/norm(dr_edge)*dr_edge;
end

function [cross_flag,stagnation_flag]=judgeAttachment(normal_vector,point_list,Ve_list)
% judge if this element have attachment line across
% 
% reference: [1] KENWRIGHT D N. Automatic detection of open and closed
% separation and attachment lines; proceedings of the Proceedings
% Visualization '98 (Cat No98CB36276), F 18-23 Oct. 1998, 1998 [C].
%
topology_torlance=1e-6;
precision_torlance=1e-12;

cross_flag=false(1);
stagnation_flag=false(0);

% project to 2D
e3=normal_vector;
e1=point_list(1,1:3)-point_list(2,1:3);
e1=e1-(e1*e3')*e3;
e1=e1/norm(e1);
e2=cross(e3,e1);

% transform local cooridinate to 2D cooridinate
% all vector just simply project to element panel
Ve_list=Ve_list*[e1;e2;e3]';
point_2D_list=point_list*[e1;e2;e3]';
point_2D_list(:,end)=1;

% interpolation fit flow
coefficient=(point_2D_list\Ve_list(:,[1,2]))';

jacobian=coefficient(:,1:2); % [b1,c1;b2,c2]
jacobian(abs(jacobian) < precision_torlance)=0; % solve precision
bias=coefficient(:,3); % [a1;a2]

% jacobian matrix is negative, process stop
det_jacobian=(jacobian(1,1)*jacobian(2,2)-jacobian(1,2)*jacobian(2,1));
if (((jacobian(1,1)+jacobian(2,2))^2-...
        4*det_jacobian) < -topology_torlance)
    eig_value=[];
    stagnation_point=[];
    return;
end

% evaluate eigenvalues of jacobian matrix
% if one of eigenvalues is zero, process stop
eig_value=eig(jacobian);
if sum(abs(eig_value) <= topology_torlance)
    stagnation_point=[];
    return;
end

stagnation_point=-jacobian\bias;
% stagnation cannot to far 
% if max(abs(stagnation_point)) > 10
%     return;
% end

% project element to canonical coordinates
[eig_vector,eig_value]=eig(jacobian);
eig_value=[eig_value(1,1),eig_value(2,2)];
if (eig_value(1) < eig_value(2))
    eig_value=fliplr(eig_value);
    eig_vector=fliplr(eig_vector);
end

% normalize point_2D
point_proj_list=(point_2D_list(:,1:2)-stagnation_point')/eig_vector';
if det(eig_vector) < 0
    point_proj_list=flipud(point_proj_list);
end

max_bou=max(abs(point_proj_list));
if ((eig_value(2)*max_bou(2))^2-(eig_value(1)*max_bou(1))^2) > cos(pi/3) % velocity do not too close
    return;
end

offset = norm(max(point_proj_list)-min(point_proj_list))*1e-2;
point_proj_list=geoCurveOffset(point_proj_list,offset);

% if (0,0) inside element
if judgeOriginSurround(point_proj_list)
    cross_flag=true(1);
    stagnation_flag=1;
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
    if judgeCrossY(point_proj_list(1:end,[1,2]))
        cross_flag=true(1);
    end
elseif ((eig_value(1) < -topology_torlance) && ...
        (eig_value(2) < -topology_torlance))
    % attracting node, check
    % concentrate line is not stagnation point
    
elseif ((eig_value(1) > topology_torlance) && (eig_value(2) < -topology_torlance))
    % saddle, judge axis X and axis Y which is separation line
    % means which eigenvalue is large than zero
    % repelling node, check small eigenvalue corresponded axis
    if judgeCrossY(point_proj_list(1:end,[1,2]))
        cross_flag=true(1);
    end
end

end

function setElementNoSpecify(element_list)
for element_index=1:length(element_list)
    if isempty(element_list(element_index).stagnation)
        element_list(element_index).stagnation=false(1);
    elseif ~element_list(element_index).stagnation
        element_list(element_index).stagnation=false(1);
    end
end
end

function setElementSpecify(element_list)
for element_index=1:length(element_list)
    element_list(element_index).stagnation=true(1);
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

function point_list=getVertexArou()
% vertex_list=user_model.vertex_list;
% vertex_empty=user_model.vertex_empty;
% 
% % identify stagnation vertex
% % base on reference [1] to extract
% % analysis element is nearby element center of vertex
% for vertex_index=1:length(vertex_list)
%     % base on around vertex to determine if stagnation vertex
%     vertex=vertex_list(vertex_index);
%     if (vertex_list(vertex_index) == vertex_empty)
%         continue;
%     end
% 
%     vertex_next_list=vertex.point_next_list;
%     element_ref_list=vertex.element_ref_list;
%     nearby_number=vertex.nearby_number;
% 
%     % average for nearby element normal vector
%     normal_vector=zeros(1,3);
%     for nearby_index=1:nearby_number
%         normal_vector=normal_vector+...
%             element_ref_list(nearby_index).normal_vector;
%     end
%     normal_vector=normal_vector/double(nearby_number);
% 
% %     if dot(normal_vector,free_flow_vector) < - 0.98
% %         disp('');
% %     end
%     
%     % only consider upwind element have surface flow
%     if (normal_vector*free_flow_vector > 0)
%         setElementNoSpecify(element_ref_list);
%         continue;
%     end
% 
%     % load data and process symmetry, make sure point index list is curve
%     % search point_index_list from element to make sure point_index_list is
%     % anti-clock wise, which make point_index_list become an visual element
%     %
%     % vertex_prev_list means: point_back->vertex_index
%     % so it can be:
%     % vertex_next_list↑
%     % vertex_index↑ element_ref_list
%     % vertex_prev_list↑
%     %
%     vertex_prev_list=zeros(nearby_number,1);
%     for nearby_index=1:nearby_number
%         point_index_list=element_ref_list(nearby_index).point_index_list;
% 
%         % find curruent point index in nearby element point index list
%         [~,index]=judgeMatExistNum(point_index_list,vertex_index);
%         index=index-1;
%         if index == 0
%             index=length(point_index_list);
%         end
%         vertex_prev_list(nearby_index)=point_index_list(index);
%     end
%     vertex.point_prev_list=vertex_prev_list;
%     
%     % base on vertex_next_list and vertex_prev_list generate vertex_arou_list
%     % for point on symmetry face, vertex_next_list and vertex_prev_list
%     % will not connect, need to process respectively
%     if ~isempty(user_model.SYMMETRY) && (abs(point_list(vertex_index,identify_dimension)) < geometry_torlance)
%         normal_vector(identify_dimension)=0;
%         normal_vector=normal_vector/norm(normal_vector);
% 
%         % process symmetry
%         % start search from first of vertex_next_list
%         % forward search
%         index_arou_list=1;
% %         vertex_arou_list=vertex_next_list(1);
%         [exit_flag,index]=judgeMatExistNum...
%             (vertex_next_list,vertex_prev_list(index_arou_list));
%         while exit_flag % if exist, add into point list
%             index_arou_list=[index_arou_list;index];
% %             vertex_arou_list=[vertex_arou_list;vertex_next_list(index_arou_list(end))];
%             [exit_flag,index]=judgeMatExistNum...
%                 (vertex_next_list,vertex_prev_list(index_arou_list(end)));
%         end
%         % if not exist, add into point list
%         index_arou_list=[index_arou_list;index];
% %         vertex_arou_list=[vertex_arou_list,vertex_prev_list(index_arou_list(end))];
% 
%         % backward search
%         [~,index]=judgeMatExistNum...
%             (vertex_prev_list,vertex_next_list(1));
%         index_arou_list=[index;index_arou_list];
%         [exit_flag,index]=judgeMatExistNum...
%             (vertex_prev_list,vertex_next_list(index_arou_list(1)));
%         while exit_flag
%             index_arou_list=[index;index_arou_list];
% %             vertex_arou_list=[vertex_prev_list(index_arou_list(1));vertex_arou_list];
%             [exit_flag,index]=judgeMatExistNum...
%                 (vertex_prev_list,vertex_next_list(index_arou_list(1)));
%         end
%         % if not exist, add into point list
%         index_arou_list=[index;index_arou_list];
% %         vertex_arou_list=[vertex_next_list(index_arou_list(1));vertex_arou_list];
% 
%         point_arou_list=zeros(length(index_arou_list),3);
%         Ve_list=zeros(length(index_arou_list),3);
%         for index=1:length(index_arou_list)
%             point_arou_list(index,:)=element_ref_list(index_arou_list(index)).center_point;
%             Ve_list(index,:)=element_ref_list(index_arou_list(index)).surface_flow;
%         end
%         
%         % because of symmetry, 2:end-1 point data need to add symmetry
%         flipur_list=flipud(Ve_list);
%         flipur_list(:,identify_dimension)=-flipur_list(:,identify_dimension);
%         Ve_list=[Ve_list;flipur_list];
%         flipur_list=flipud(point_arou_list);
%         flipur_list(:,identify_dimension)=-flipur_list(:,identify_dimension);
%         point_arou_list=[point_arou_list;flipur_list];
%     else
%         % start search from first of vertex_next_list
%         index_arou_list=zeros(nearby_number,1);
%         index_arou_list(1)=1;
% %         vertex_arou_list=zeros(nearby_number,1);
% %         vertex_arou_list(1)=vertex_next_list(index_arou_list(1));
%         for nearby_index=2:nearby_number
%             [~,index_arou_list(nearby_index)]=judgeMatExistNum...
%                 (vertex_next_list,vertex_prev_list(index_arou_list(nearby_index-1)));
% %             vertex_arou_list(nearby_index)=vertex_next_list(index_arou_list(nearby_index));
%         end
%         
%         point_arou_list=zeros(length(index_arou_list),3);
%         Ve_list=zeros(length(index_arou_list),3);
%         for index=1:length(index_arou_list)
%             point_arou_list(index,:)=element_ref_list(index_arou_list(index)).center_point;
%             Ve_list(index,:)=element_ref_list(index_arou_list(index)).surface_flow;
%         end
%     end
% 
%     [cross_flag,eig_value,stagnation_point]=judgeSpecify(normal_vector,point_arou_list,Ve_list);
%     if cross_flag
%         setElementSpecify(element_ref_list)
%         fprintf('eig: %f,%f,  stag: %f,%f\n',eig_value,stagnation_point)
%     else
%         setElementNoSpecify(element_ref_list)
%     end
% end
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
