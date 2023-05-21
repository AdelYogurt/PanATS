function [part_list,point_list,geometry]=convertSTLToMesh...
    (part_list,remove_redundance,geometry_torlance)
% convert STL format part into mesh format part
%
% if remove_redundance, function will call ADTreePoint to remove repeat
% point, change element_index to same coordinate point
%
if nargin < 2
    remove_redundance=true(1);
end

if ~iscell(part_list)
    part_list={part_list};
end

if nargin < 3 || isempty(geometry_torlance)
    geometry_torlance=1e-12;
end

% add all stl element_list to point list
point_list=[];
for part_index=1:length(part_list)
    mesh_list=part_list{part_index}.mesh_list;

    % process each mesh
    for mesh_index=1:length(mesh_list)
        mesh=mesh_list{mesh_index};
        if ~strcmp(mesh.element_type,'stl')
            error('convertSTLToMesh: element_type of mesh is not stl format');
        end

        % add point element_list into point_list
        point_list=[point_list;mesh.element_list];
    end
end

if remove_redundance
    % delete same coordination point
    % need ADTreePoint
    [ADT,index_list]=ADTreePoint(point_list,[],[],geometry_torlance);

    exist_point_number=0; % offset index of point index list
    for part_index=1:length(part_list)
        part=part_list{part_index};
        mesh_list=part.mesh_list;

        % process each mesh
        for mesh_index=1:length(mesh_list)
            mesh=mesh_list{mesh_index};
            element_list=mesh.element_list;

            [mesh_point_number,~]=size(element_list);
            mesh_point_number=int32(mesh_point_number);
            element_number=mesh_point_number/3;

            % convert element_list
            element_list=index_list(exist_point_number+1:exist_point_number+mesh_point_number);
            element_list=reshape(int32(element_list),3,element_number)';
            
            % sort element
            mesh.element_list=element_list;
            mesh.element_type='S3';
            mesh.element_ID=int8(5);

            mesh_list{mesh_index}=mesh;

            exist_point_number=exist_point_number+mesh_point_number;
        end

        part.mesh_list=mesh_list;
        part_list{part_index}=part;
    end

    geometry.min_bou=ADT.min_bou;
    geometry.max_bou=ADT.max_bou;
    geometry.dimension=3;
else
    % do not delete same coordination point
    exist_point_number=0; % offset index of point index list
    for part_index=1:length(part_list)
        part=part_list{part_index};
        mesh_list=part.mesh_list;

        % process each mesh
        for mesh_index=1:length(mesh_list)
            mesh=mesh_list{mesh_index};
            if ~strcmp(mesh.element_type,'stl')
                error('convertSTLToMesh: element_type of mesh is not stl format');
            end

            element_list=mesh.element_list;

            [mesh_point_number,~]=size(element_list);
            mesh_point_number=int32(mesh_point_number);

            % convert element list
            element_number=mesh_point_number/3;
            element_list=reshape(int32(1:mesh_point_number),3,element_number)'+exist_point_number;
            
            % sort element
            mesh.element_list=element_list;
            mesh.element_type='S3';
            mesh.element_ID=int8(5);

            mesh_list{mesh_index}=mesh;

            exist_point_number=exist_point_number+mesh_point_number;
        end

        part.mesh_list=mesh_list;
        part_list{part_index}=part;
    end

    geometry.min_bou=min(point_list);
    geometry.max_bou=max(point_list);
    geometry.dimension=3;
end

if length(part_list) == 1
    part_list=part_list{1};
end

end
