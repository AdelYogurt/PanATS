function mesh_data=convertWGSToMesh(mesh_data,remove_redundance,geometry_torlance)
% convert LaWGS format part into mesh format part
%
% if remove_redundance, function will call ADTreePoint to remove repeat
% point, change element_idx to same coordinate point
%
% input:
% mesh_data, remove_redundance, geometry_torlance
%
% output:
% mesh_data
%
if nargin < 2 || isempty(remove_redundance)
    remove_redundance=true(1);
end
if nargin < 3 || isempty(geometry_torlance)
    geometry_torlance=1e-12;
end

marker_name_list=fieldnames(mesh_data);

% add all surface to point list
point_list=[];
for marker_index=1:length(marker_name_list)
    marker_name=marker_name_list{marker_index};
    if strcmp(marker_name,'geometry')
        continue;
    end

    marker=mesh_data.(marker_name);
    if ~strcmp(marker.type,'wgs')
        error('convertWGSToSTL: type of mesh is not LaWGS format');
    end

    X=marker.X;
    Y=marker.Y;
    Z=marker.Z;

    if any(size(X) ~= size(Y)) || any(size(X) ~= size(Z))
        error('convertWGSToSTL: size of X,Y,Z in mesh are not equal');
    end

    % add this mesh point into point list
    point_list=[point_list;[X(:),Y(:),Z(:)]];
end

if remove_redundance
    % delete same coordination point
    % need ADTreePoint
    [~,index_list]=ADTreePoint(point_list,[],[],geometry_torlance);
    [index_list,~,map_list]=unique(index_list);
    point_list=point_list(index_list,:);

    node_idx=uint32(0); % offset idx of point idx list
    for marker_index=1:length(marker_name_list)
        marker_name=marker_name_list{marker_index};
        if strcmp(marker_name,'geometry')
            continue;
        end
        marker=mesh_data.(marker_name);

        [point_number,line_number]=size(marker.X);

        % initialize element sort space
        calculate_element=(point_number-1)*(line_number-1)*2;
        element_indx=0;
        element_list=zeros(calculate_element,3);

        % sort element
        for line_idx=1:line_number-1
            base_offset=point_number*(line_idx-1);
            for point_idx=1:point_number-1
                point1_idx=map_list(base_offset+point_idx+node_idx);
                point2_idx=map_list(base_offset+point_idx+point_number+node_idx);
                point3_idx=map_list(base_offset+point_idx+point_number+node_idx+1);
                point4_idx=map_list(base_offset+point_idx+node_idx+1);

                if point1_idx == point2_idx || point2_idx == point4_idx || point4_idx == point1_idx
                    % small element degeneration to line, discard it
                else
                    element_indx=element_indx+1;
                    element_list(element_indx,:)=[uint32(point1_idx),uint32(point2_idx),uint32(point4_idx)];
                end

                if point2_idx == point3_idx || point3_idx == point4_idx || point4_idx == point2_idx
                    % small element degeneration to line, discard it
                else
                    element_indx=element_indx+1;
                    element_list(element_indx,:)=[uint32(point2_idx),uint32(point3_idx),uint32(point4_idx)];
                end
            end
        end
        node_idx=node_idx+uint32(point_number*line_number);

        marker=rmfield(marker,'X');
        marker=rmfield(marker,'Y');
        marker=rmfield(marker,'Z');

        % sort element
        marker.type='TRI_3';
        marker.ID=uint8(5);
        marker.element_list=element_list(1:element_indx,:);
        marker.number_list=3;

        mesh_data.(marker_name)=marker;
    end
else
    % do not delete same coordination point
    node_idx=uint32(0); % offset idx of point idx list

    for marker_index=1:length(marker_name_list)
        marker_name=marker_name_list{marker_index};
        if strcmp(marker_name,'geometry')
            continue;
        end
        marker=mesh_data.(marker_name);

        [point_number,line_number]=size(marker.X);

        % initialize element sort space
        calculate_element=(point_number-1)*(line_number-1)*2;
        element_indx=0;
        element_list=zeros(calculate_element,3);

        % sort element
        for line_idx=1:line_number-1
            base_offset=point_number*(line_idx-1);
            for point_idx=1:point_number-1
                point1_idx=base_offset+point_idx+node_idx;
                point2_idx=base_offset+point_idx+point_number+node_idx;
                point3_idx=base_offset+point_idx+point_number+node_idx+1;
                point4_idx=base_offset+point_idx+node_idx+1;

                d12=point_list(point2_idx,:)-point_list(point1_idx,:);
                d24=point_list(point4_idx,:)-point_list(point2_idx,:);

                if norm(d12) < eps || norm(d24) < eps
                    % small element degeneration to line, discard it
                else
                    element_indx=element_indx+1;
                    element_list(element_indx,:)=[uint32(point1_idx),uint32(point2_idx),uint32(point4_idx)];
                end

                d23=point_list(point3_idx,:)-point_list(point2_idx,:);
                d34=point_list(point4_idx,:)-point_list(point3_idx,:);

                if norm(d23) < eps || norm(d34) < eps
                    % small element degeneration to line, discard it
                else
                    element_indx=element_indx+1;
                    element_list(element_indx,:)=[uint32(point2_idx),uint32(point3_idx),uint32(point4_idx)];
                end
            end
        end
        node_idx=node_idx+uint32(point_number*line_number);

        marker=rmfield(marker,'X');
        marker=rmfield(marker,'Y');
        marker=rmfield(marker,'Z');

        % sort element
        marker.type='TRI_3';
        marker.ID=uint8(5);
        marker.element_list=element_list(1:element_indx,:);
        marker.number_list=3;

        mesh_data.(marker_name)=marker;
    end
end

geometry.point_list=point_list;
geometry.dimension=3;
mesh_data.geometry=geometry;

end
