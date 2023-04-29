function writeMeshWGS(filename_mesh,part_list)
% write LaWGS format mesh data into wgs file
%
% input:
% filename_mesh, part_list 
%
% notice:
% point_list is coordinate of all node
% part_list{part.name, part.mesh_list{mesh.element_list, mesh.element_type, mesh.element_ID}}
%
if strcmp(filename_mesh(end-3:end),'.wgs')
    filename_mesh=filename_mesh(1:end-4);
end

file_mesh=fopen([filename_mesh,'.wgs'],'w');

if ~iscell(part_list)
    part_list={part_list};
end

% write basic information of inp file
fprintf(file_mesh,'Create by writeMeshWGS.m\n');

% write start of part
for part_index=1:length(part_list)
    part=part_list{part_index};
    mesh_list=part.mesh_list;

    % write part name
    if length(mesh_list) > 1
        fprintf(file_mesh,'%s1\n',part.name);
    else
        fprintf(file_mesh,'%s\n',part.name);
    end

    mesh=mesh_list{1};
    X=mesh.X;Y=mesh.Y;Z=mesh.Z;
    [point_number,line_number]=size(X);

    % write part information data
    fprintf(file_mesh,' %d %d %d 0 0 0 0 0 0 0 1 1 1 0\n',part_index,line_number,point_number);

    % write element data
    for line_index=1:line_number
        for point_index=1:point_number
            fprintf(file_mesh,'%14f %14f %14f\n',X(point_index,line_index),Y(point_index,line_index),Z(point_index,line_index));
        end
    end

    if length(mesh_list) > 1
        for mesh_index=2:length(mesh_list)
            % write part name
            fprintf(file_mesh,'%s%d\n',part.name,mesh_index);

            mesh=mesh_list{mesh_index};
            X=mesh.X;Y=mesh.Y;Z=mesh.Z;
            [point_number,line_number]=size(X);

            % write part information data
            fprintf(file_mesh,' %d %d %d 0 0 0 0 0 0 0 1 1 1 0\n',part_index,line_number,point_number);

            % write element data
            for line_index=1:line_number
                for point_index=1:point_number
                    fprintf(file_mesh,'%14f %14f %14f\n',X(point_index,line_index),Y(point_index,line_index),Z(point_index,line_index));
                end
            end
        end
    end
end

fclose(file_mesh);
clear('file_mesh');

end