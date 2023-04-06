function [point_list,part_list]=convertWGSToMesh(part_list,remove_redundance)
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

if remove_redundance

else
    % do not delete same coordination point
    point_list=[];
    for part_index=1:length(part_list)
        part=part_list{part_index};
        mesh_list=part.mesh_list;

        % process each mesh
        for mesh_index=1:length(mesh_list)
            exist_point_number=size(point_list,1);

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

            [point_number,line_number]=size(X);
            point_mesh_list=zeros(point_number*line_number,3);

            % add point of X, Y, Z into point_list
            for line_index=1:line_number
                base_offset=point_number*(line_index-1);
                for point_index=1:point_number
                    point_mesh_list(base_offset+point_index,:)=...
                        [X(point_index,line_index),Y(point_index,line_index),Z(point_index,line_index)];
                end
            end

            % write each element
            calculate_element=(point_number-1)*(line_number-1)*2;
            element_number=0;
            element_list=zeros(calculate_element,3);
            
            % write element
            for line_index=1:line_number-1
                base_offset=point_number*(line_index-1);
                for point_index=1:point_number-1
                    point1_index=base_offset+point_index;
                    point2_index=base_offset+point_index+1;
                    point3_index=base_offset+point_index+1+point_number;
                    point4_index=base_offset+point_index+point_number;

                    d12=point_mesh_list(point2_index,:)-point_mesh_list(point1_index,:);
                    d23=point_mesh_list(point3_index,:)-point_mesh_list(point2_index,:);

                    if norm(d12) < eps || norm(d23) < eps
                        % small element degeneration to line, discard it
                    else
                        element_number=element_number+1;
                        element_list(element_number,:)=...
                            [int32(point1_index+exist_point_number),int32(point2_index+exist_point_number),int32(point3_index+exist_point_number)];
                    end

                    d34=point_mesh_list(point4_index,:)-point_mesh_list(point3_index,:);
                    d41=point_mesh_list(point1_index,:)-point_mesh_list(point4_index,:);

                    if norm(d34) < eps || norm(d41) < eps
                        % small element degeneration to line, discard it
                    else
                        element_number=element_number+1;
                        element_list(element_number,:)=...
                            [int32(point3_index+exist_point_number),int32(point4_index+exist_point_number),int32(point1_index+exist_point_number)];
                    end
                end
            end

            element_list=element_list(1:element_number,:);
            element_type='S3';
            element_ID=int8(5);

            mesh=rmfield(mesh,'X');
            mesh=rmfield(mesh,'Y');
            mesh=rmfield(mesh,'Z');

            mesh.element_list=element_list;
            mesh.element_type=element_type;
            mesh.element_ID=element_ID;

            mesh_list{mesh_index}=mesh;

            % add this mesh point into point list
            point_list=[point_list;point_mesh_list];
        end

        part.mesh_list=mesh_list;
        part_list{part_index}=part;
    end
end

end