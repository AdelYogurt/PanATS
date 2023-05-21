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

        dr_norm=norm(dr);
        dr_ref_norm=norm(dr_ref);

        if dr_ref_norm == 0 || dr_norm == 0
            continue;
        end

        % check if edge was upstream 
        if ((dr/dr_norm+dr_ref/dr_ref_norm)*surface_flow') <= -1e-12 % float error
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
            if ((dedge_ref/norm(dedge_ref)+dedge_prev/norm(dedge_prev))*element_ref.surface_flow') <= -1e-12 % float error
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
        % search next edge
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

                if isempty(point_end)
                    point_end=element.center_point;
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
        if ((dr/norm(dr)+dr_ref/norm(dr_ref))*surface_flow') >= 1e-12 % float error
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
        if ((dr/norm(dr)+dr_ref/norm(dr_ref))*surface_flow') >= 1e-12 % float error
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
        element.streamline_ref=[element.streamline_ref;...
            streamline_index,index];
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
            if ((dedge_ref/norm(dedge_ref)+dedge_prev/norm(dedge_prev))*element_ref.surface_flow') >= 1e-12 % float error
                % check if edge was crossed by surface_flow
                if cross(element_ref.surface_flow,dedge_ref)*cross(element_ref.surface_flow,dedge_prev)' <= 0
                    element_search_flag=true(1);
                    break; % goto edge search
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

    if isempty(element_ref) || sum(abs(element_ref.surface_flow)) < geometry_torlance
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
