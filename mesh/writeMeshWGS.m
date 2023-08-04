function writeMeshWGS(mesh_filestr,grid,marker_name_list)
% write LaWGS format mesh data into wgs file
%
% input:
% filename_mesh, part_list 
%
% notice:
% point_list is coordinate of all node
% part_list{part.name, part.mesh_list{mesh.element_list, mesh.element_type, mesh.element_ID}}
%
if nargin < 3
    marker_name_list=[];
end
if isempty(marker_name_list)
    marker_name_list=fieldnames(grid);
end

% check file name
if length(mesh_filestr) > 4
    if ~strcmpi(mesh_filestr((end-3):end),'.wgs')
        mesh_filestr=[mesh_filestr,'.wgs'];
    end
else
    mesh_filestr=[mesh_filestr,'.wgs'];
end

file_mesh=fopen(mesh_filestr,'w');

% write basic information of inp file
fprintf(file_mesh,'Create by writeMeshWGS.m\n');

% write each marker to file
for marker_index=1:length(marker_name_list)
    marker_name=marker_name_list{marker_index};
    if strcmp(marker_name,'point_list') || strcmp(marker_name,'name')
        continue;
    end

    % write part name
    fprintf(file_mesh,'''%s''\n',marker_name);

    X=grid.(marker_name).X;Y=grid.(marker_name).Y;Z=grid.(marker_name).Z;
    [point_number,line_number]=size(X);

    % write part information data
    fprintf(file_mesh,' %d %d %d 0 0 0 0 0 0 0 1 1 1 0\n',marker_index,line_number,point_number);

    % write element data
    for line_index=1:line_number
        for point_index=1:point_number
            fprintf(file_mesh,'%14f %14f %14f\n',X(point_index,line_index),Y(point_index,line_index),Z(point_index,line_index));
        end
    end

end

fclose(file_mesh);
clear('file_mesh');

end