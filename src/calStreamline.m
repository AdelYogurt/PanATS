function [streamline,line_len,cross_idx]=calStreamline(elem_idx,...
    point_list,element_list,boolean_attach_list,elem_flow_list,normal_vector_list)
% calculate streamline from point_start and element
%
geometry_torlance=1e-9;

elem_num=length(element_list);
max_size=elem_num;

streamline=zeros(max_size,3); % cross point
line_len=zeros(max_size,1); % current length
cross_idx=zeros(max_size,1,'uint32'); % [cross elem_idx, cross edge_idx]
data_idx=1;
point_start=sum(point_list(element_list(elem_idx).Node_idx,:),1)/element_list(elem_idx).node_num;
% line(point_start(:,1),point_start(:,2),point_start(:,3),'Marker','o','Color','b')
streamline(data_idx,:)=point_start;
cross_idx(data_idx,:)=elem_idx;

%% upwind search begin

while elem_idx ~= 0 && data_idx <= max_size/2
    % load element information
    elem=element_list(elem_idx);
    Node_idx=elem.Node_idx;
    node_num=elem.node_num;
    elem_flow=-elem_flow_list(elem_idx,:);
    normal_vector=normal_vector_list(elem_idx,:);
    norm_elem_flow=norm(elem_flow);
    if norm_elem_flow <= geometry_torlance % end in no flow element
        break;
    end
    elem_flow=elem_flow/norm_elem_flow;
    node_list=point_list(Node_idx,:);

    % serach flow upwind end point
    [point_end,node_idx,d_len]=calElementCross...
        (point_start,node_list,node_num,normal_vector,elem_flow,geometry_torlance);
    if d_len < geometry_torlance % try upwind by node point
        proj_direct=(node_list-point_start)*elem_flow';
        [proj_direct_max,node_max_idx]=max(proj_direct);node_max_idx=node_max_idx(1);
        if proj_direct_max > 0
            point_end=node_list(node_max_idx,:);
            d_len=norm(point_start-point_end);
        end
    end

    if d_len < geometry_torlance
        % cannot find upwind end point, stop
        break;
    end

    % record parameter
    data_idx=data_idx+1;
    streamline(data_idx,:)=point_end;
    line_len(data_idx)=line_len(data_idx-1)+d_len;
    cross_idx(data_idx)=elem_idx;
    % line(point_end(:,1),point_end(:,2),point_end(:,3),'Marker','o','Color','r')

    % end in attach element
    if boolean_attach_list(elem_idx)
        break;
    end
   
    % search next element
    % check if is point_start overlap node list
    node_overlap_idx=find(vecnorm(node_list-point_end,2,2) < geometry_torlance,1);
    if ~isempty(node_overlap_idx)
        % if overlap, we need to find next element that can upwind
        [Elem_arou,Vertex_arou]=getAroundElement(element_list,elem_idx,node_overlap_idx);
        Direct_arou=(point_list(Vertex_arou,:)-point_end);
        Direct_arou=Direct_arou./vecnorm(Direct_arou,2,2);
        Angle_arou=Direct_arou*elem_flow';
        [~,upwind_idx]=max(Angle_arou);
        elem_idx=Elem_arou(upwind_idx);
    else
        % if not overlap, use node_idx to find next element to upwind
        elem_idx_new=elem.Vertex_next(node_idx);
        if elem_idx_new ~= 0,elem_idx=elem_idx_new;end
    end
    point_start=point_end;
end


%% reverse and connect point

if data_idx > 1
    % reverse streamline
    streamline(1:data_idx,:)=flipud(streamline(1:data_idx,:));
    line_len(1:data_idx)=flipud(line_len(1:data_idx));
    line_len(1:data_idx)=line_len(1)-line_len(1:data_idx);

    % reverse cross_list
    cross_idx(1:data_idx)=flipud(cross_idx(1:data_idx));
    cross_idx(2:data_idx)=cross_idx(1:data_idx-1);

    data_idx=data_idx-1;
    point_start=streamline(data_idx,1:3);
    elem_idx=cross_idx(data_idx+1);
end

%% downwind search begin

while elem_idx ~= 0 && data_idx <= max_size
    % load element information
    elem=element_list(elem_idx);
    Node_idx=elem.Node_idx;
    node_num=elem.node_num;
    elem_flow=elem_flow_list(elem_idx,:);
    normal_vector=normal_vector_list(elem_idx,:);
    norm_elem_flow=norm(elem_flow);
    if norm_elem_flow <= geometry_torlance % end in no flow element
        break;
    end
    elem_flow=elem_flow/norm_elem_flow;
    node_list=point_list(Node_idx,:);

    % calculate next cross point
    [point_end,node_idx,d_len]=calElementCross...
        (point_start,node_list,node_num,normal_vector,elem_flow,geometry_torlance);

    if d_len < geometry_torlance
        % check whether point_start local near node point cause
        Bool_overlap=vecnorm(point_start-node_list,2,2) < geometry_torlance;
        if any(Bool_overlap)
            % if yes, search element that can continue search
            elem_idx_origin=elem_idx;
            node_idx=find(Bool_overlap);
            vertex_idx=Node_idx(node_idx);vertex_idx=vertex_idx(1); % overlap point

            % base on half edge to search all around element
            elem_idx=elem.Vertex_next(node_idx);
            while elem_idx_origin ~= elem_idx
                elem=element_list(elem_idx);
                Node_idx=elem.Node_idx;
                node_num=elem.node_num;
                elem_flow=elem_flow_list(elem_idx,:);
                normal_vector=normal_vector_list(elem_idx,:);
                norm_elem_flow=norm(elem_flow);
                if norm_elem_flow <= geometry_torlance % end in no flow element
                    break;
                end
                elem_flow=elem_flow/norm_elem_flow;
                node_list=point_list(Node_idx,:);

                [point_end,node_idx,d_len]=calElementCross...
                    (point_start,node_list,node_num,normal_vector,elem_flow,geometry_torlance);

                if d_len >= geometry_torlance
                    % next element found, contiune searching
                    break;
                else
                    % try next point
                    node_idx=find(Node_idx == vertex_idx);
                    elem_idx=elem.Vertex_next(node_idx);
                    if elem_idx == 0
                        break;
                    end
                end
            end

            if elem_idx_origin == elem_idx || elem_idx == 0
                % if all around element searched but there is not element found
                break;
            end
        else
            % if no, means encounter the opposite flow, stop search
            break;
        end
    end

    % record parameter
    data_idx=data_idx+1;
    streamline(data_idx,:)=point_end;
    line_len(data_idx)=line_len(data_idx-1)+d_len;
    cross_idx(data_idx)=elem_idx;
    %     line(point_end(:,1),point_end(:,2),point_end(:,3),'Marker','o','Color','g')

    % search next element
    elem_idx=elem.Vertex_next(node_idx);
    point_start=point_end;
end

streamline=streamline(1:data_idx,:);
line_len=line_len(1:data_idx,:);
cross_idx=cross_idx(1:data_idx);
end

function [point_end,node_idx,d_len]=calElementCross...
    (point_start,node_list,node_num,normal_vector,elem_flow,geometry_torlance)
% calculate next cross point
%

% elem_flow as x, normal_vector as z, project point
y_flow=cross(normal_vector,elem_flow);
rotation=[elem_flow;y_flow;normal_vector]';
proj_node_list=(node_list-point_start)*rotation;

% base on cross x axis to calculate projected point_end
point_end=[min(proj_node_list(:,1)),0,0];
node_idx=0;
for node_base_idx=1:node_num
    if node_base_idx == node_num
        node_ref_idx=1;
    else
        node_ref_idx=node_base_idx+1;
    end

    vertex_base_x=proj_node_list(node_base_idx,1);
    vertex_ref_x=proj_node_list(node_ref_idx,1);
    vertex_base_y=proj_node_list(node_base_idx,2);
    vertex_ref_y=proj_node_list(node_ref_idx,2);

    if (vertex_base_y <= geometry_torlance && vertex_ref_y >= -geometry_torlance) ||...
            (vertex_base_y >= -geometry_torlance && vertex_ref_y <= geometry_torlance)
        if vertex_base_y == vertex_ref_y
            if vertex_ref_x > point_end(1)
                point_end=[vertex_ref_x,0,0];
                node_idx=node_base_idx;
            end
        else
            cross_point=vertex_base_x-vertex_base_y*(vertex_base_x-vertex_ref_x)/(vertex_base_y-vertex_ref_y);
            if cross_point > point_end(1)
                point_end=[cross_point,0,0];
                node_idx=node_base_idx;
            end
        end
    end
end

if node_idx
    % re project to real coordinate
    point_end=point_end*rotation'+point_start;
    d_len=norm(point_start-point_end);
else
    d_len=0;
    point_end=[];
end

end

