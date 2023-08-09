function mesh_data=convertSTLToMesh(mesh_data,remove_redundance,geometry_torlance)
% convert STL format part into mesh format part
%
% if remove_redundance, function will call ADTreePoint to remove repeat
% point, change element_index to same coordinate point
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
    if ~strcmp(marker.type,'stl')
        error('convertSTLToMesh: element_type of mesh is not stl format');
    end

    % add point element_list into point_list
    point_list=[point_list;marker.element_list];
end

if remove_redundance
    % delete same coordination point
    % need ADTreePoint
    [~,index_list]=ADTreePoint(point_list,[],[],geometry_torlance);
    [index_list,~,map_list]=unique(index_list);
    point_list=point_list(index_list,:);

    node_index=uint32(0); % offset index of point index list
    for marker_index=1:length(marker_name_list)
        marker_name=marker_name_list{marker_index};
        if strcmp(marker_name,'geometry')
            continue;
        end
        marker=mesh_data.(marker_name);

        % process each mesh
        element_list=marker.element_list;

        marker_point_number=size(element_list,1);
        marker_point_number=uint32(marker_point_number);
        marker_element_number=marker_point_number/3;

        % convert element_list
        element_list=map_list(node_index+1:node_index+marker_point_number);
        element_list=reshape(uint32(element_list),3,marker_element_number)';

        node_index=node_index+uint32(marker_point_number);

        % sort element
        marker.type='TRI_3';
        marker.ID=uint8(5);
        marker.element_list=element_list;
        mesh_data.(marker_name).number_list=3;

        mesh_data.(marker_name)=marker;
    end

else
    % do not delete same coordination point
    node_index=uint32(0); % offset idx of point idx list

    for marker_index=1:length(marker_name_list)
        marker_name=marker_name_list{marker_index};
        if strcmp(marker_name,'geometry')
            continue;
        end
        marker=mesh_data.(marker_name);

        element_list=mesh_data.(marker_name).element_list;

        marker_point_number=size(element_list,1);
        marker_point_number=uint32(marker_point_number);
        marker_element_number=marker_point_number/3;

        element_list=reshape(uint32(1:marker_point_number),3,marker_element_number)'+node_index;
        node_index=node_index+marker_point_number;

        % sort element
        marker.type='TRI_3';
        marker.ID=uint8(5);
        marker.element_list=element_list;
        marker.number_list=3;

        mesh_data.(marker_name)=marker;
    end
end

geometry.point_list=point_list;
geometry.dimension=3;
mesh_data.geometry=geometry;

end
