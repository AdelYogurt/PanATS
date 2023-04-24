function [part_list,point_list,geometry]=convertWGSToMesh...
    (part_list,remove_redundance,geometry_torlance)
% convert LaWGS format part into mesh format part
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

% add all surface to point list
point_list=[];
for part_index=1:length(part_list)
    mesh_list=part_list{part_index}.mesh_list;

    % process each mesh
    for mesh_index=1:length(mesh_list)
        mesh=mesh_list{mesh_index};
        if ~strcmp(mesh.element_type,'wgs')
            error('convertWGSToSTL: element_type of mesh is not LaWGS format');
        end

        X=mesh.X;
        Y=mesh.Y;
        Z=mesh.Z;

        if any(size(X) ~= size(Y)) || any(size(X) ~= size(Z))
            error('convertWGSToSTL: size of X,Y,Z in mesh are not equal');
        end

        point_mesh_list=[X(:),Y(:),Z(:)];

        % add this mesh point into point list
        point_list=[point_list;point_mesh_list];
    end
end

if remove_redundance
    % delete same coordination point
    % need ADTreePoint
    [ADT,index_list]=ADTreePoint(point_list,[],[],geometry_torlance);
    
    exist_point_number=int32(0); % offset index of point index list
    for part_index=1:length(part_list)
        mesh_list=part_list{part_index}.mesh_list;
        part_element_number=0;
        
        % process each mesh
        for mesh_index=1:length(mesh_list)
            mesh=mesh_list{mesh_index};

            [point_number,line_number]=size(mesh.X);
            mesh_point_number=point_number*line_number;

            % initialize element sort space
            calculate_element=(point_number-1)*(line_number-1)*2;
            element_number=0;
            element_list=zeros(calculate_element,3);
            
            % sort element(S3)
            for line_index=1:line_number-1
                base_offset=point_number*(line_index-1);
                for point_index=1:point_number-1
                    point1_index=index_list(base_offset+point_index+exist_point_number);
                    point2_index=index_list(base_offset+point_index+1+exist_point_number);
                    point3_index=index_list(base_offset+point_index+1+point_number+exist_point_number);
                    point4_index=index_list(base_offset+point_index+point_number+exist_point_number);

                    if point1_index == point2_index || point2_index == point3_index || point1_index == point3_index
                        % small element degeneration to line, discard it
                    else
                        element_number=element_number+1;
                        element_list(element_number,:)=...
                            [int32(point1_index),int32(point2_index),int32(point3_index)];
                    end

                    if point3_index == point4_index || point4_index == point1_index || point1_index == point3_index
                        % small element degeneration to line, discard it
                    else
                        element_number=element_number+1;
                        element_list(element_number,:)=...
                            [int32(point3_index),int32(point4_index),int32(point1_index)];
                    end
                end
            end

            mesh=rmfield(mesh,'X');
            mesh=rmfield(mesh,'Y');
            mesh=rmfield(mesh,'Z');

            % sort element
            mesh.element_list=element_list(1:element_number,:);
            mesh.element_type='S3';
            mesh.element_ID=int8(5);

            mesh_list{mesh_index}=mesh;

            exist_point_number=exist_point_number+int32(mesh_point_number);
            part_element_number=part_element_number+element_number;
        end

        part_list{part_index}.mesh_list=mesh_list;
        part_list{part_index}.element_number=part_element_number;
    end

    geometry.min_bou=ADT.min_bou;
    geometry.max_bou=ADT.max_bou;
    geometry.dimension=3;
else
    % do not delete same coordination point
    exist_point_number=int32(0); % offset index of point index list
    for part_index=1:length(part_list)
        mesh_list=part_list{part_index}.mesh_list;

        % process each mesh
        for mesh_index=1:length(mesh_list)
            mesh=mesh_list{mesh_index};

            [point_number,line_number]=size(mesh.X);
            mesh_point_number=point_number*line_number;

            % initialize element sort space
            calculate_element=(point_number-1)*(line_number-1)*2;
            element_number=0;
            element_list=zeros(calculate_element,3);
            
            % sort element(S3)
            for line_index=1:line_number-1
                base_offset=point_number*(line_index-1);
                for point_index=1:point_number-1
                    point1_index=base_offset+point_index+exist_point_number;
                    point2_index=base_offset+point_index+1+exist_point_number;
                    point3_index=base_offset+point_index+1+point_number+exist_point_number;
                    point4_index=base_offset+point_index+point_number+exist_point_number;

                    d12=point_list(point2_index,:)-point_list(point1_index,:);
                    d23=point_list(point3_index,:)-point_list(point2_index,:);

                    if norm(d12) < eps || norm(d23) < eps
                        % small element degeneration to line, discard it
                    else
                        element_number=element_number+1;
                        element_list(element_number,:)=...
                            [int32(point1_index),int32(point2_index),int32(point3_index)];
                    end

                    d34=point_list(point4_index,:)-point_list(point3_index,:);
                    d41=point_list(point1_index,:)-point_list(point4_index,:);

                    if norm(d34) < eps || norm(d41) < eps
                        % small element degeneration to line, discard it
                    else
                        element_number=element_number+1;
                        element_list(element_number,:)=...
                            [int32(point3_index),int32(point4_index),int32(point1_index)];
                    end
                end
            end

            mesh=rmfield(mesh,'X');
            mesh=rmfield(mesh,'Y');
            mesh=rmfield(mesh,'Z');

            % sort element
            mesh.element_list=element_list(1:element_number,:);
            mesh.element_type='S3';
            mesh.element_ID=int8(5);

            mesh_list{mesh_index}=mesh;

            exist_point_number=exist_point_number+int32(mesh_point_number);
        end

        part_list{part_index}.mesh_list=mesh_list;
    end

    geometry.min_bou=min(point_list);
    geometry.max_bou=max(point_list);
    geometry.dimension=3;
end

if length(part_list) == 1
    part_list=part_list{1};
end

end
