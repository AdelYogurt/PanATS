function writeMeshObj(mesh_data,mesh_filestr,marker_name_list)
% write Binary STL mesh file
%
% input:
% mesh_filestr, mesh_data, marker_name_list(default all markers)
%
% notice:
% mesh.element_ID if not is 5(tri), will be convert into tri
% part(part.name, part.mesh_list{mesh.element_list, mesh.element_type, mesh.element_ID})
%
if nargin < 3
    marker_name_list=[];
    if nargin < 2
        mesh_filestr=[];
    end
end
if isempty(mesh_filestr),mesh_filestr='mesh.obj';end
if isempty(marker_name_list),marker_name_list=fieldnames(mesh_data);end
marker_index=1;
while marker_index <= length(marker_name_list)
    marker_name=marker_name_list{marker_index};
    if strcmp(marker_name,'geometry'),marker_name_list(marker_index)=[];
    else,marker_index=marker_index+1;end
end

mesh_file=fopen(mesh_filestr,'w');

dimension=mesh_data.geometry.dimension;

% write point_list
point_list=mesh_data.geometry.point_list;
for point_index=1:size(point_list,1)
    fprintf(mesh_file,'v');
    fprintf(mesh_file,' %e',point_list(point_index,1:dimension));
    fprintf(mesh_file,'\n');
end

% write element_list
for marker_index=1:length(marker_name_list)
    marker_name=marker_name_list{marker_index};
    marker=mesh_data.(marker_name);

    fprintf(mesh_file,'g %s\n',marker_name);
    ID=marker.ID;
    element_list=marker.element_list;
    if ID == 20 % mixed element
        data_index=1;
        for element_index=1:length(number_list)
            node_num=number_list(element_index);
            fprintf(mesh_file,'f');
            fprintf(mesh_file,' %d',element_list((data_index+1):(data_index+node_num)));
            fprintf(mesh_file,'\n');
            data_index=data_index+node_num+1;
        end
    else
        for element_index=1:size(element_list,1)
            fprintf(mesh_file,'f');
            fprintf(mesh_file,' %d',element_list(element_index,:));
            fprintf(mesh_file,'\n');
        end
    end
end

fclose(mesh_file);
clear('mesh_file');

end