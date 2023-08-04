function writeMeshSTL(mesh_filestr,grid,marker_name_list)
% write Binary STL mesh file
%
% input:
% filename_mesh, part(only one part)
%
% notice:
% mesh.element_ID if not is 5(tri), will be convert into tri
% part(part.name, part.mesh_list{mesh.element_list, mesh.element_type, mesh.element_ID})
%
if nargin < 3
    marker_name_list=[];
end
if isempty(marker_name_list)
    marker_name_list=fieldnames(grid);
end

% check file name
if length(mesh_filestr) > 4
    if ~strcmpi(mesh_filestr((end-3):end),'.stl')
        mesh_filestr=[mesh_filestr,'.stl'];
    end
else
    mesh_filestr=[mesh_filestr,'.stl'];
end

mesh_file=fopen(mesh_filestr,'w');

% load mesh name
mesh_name=grid.name;

% write name to file
name_length=length(mesh_name);
fwrite(mesh_file,['solid ',mesh_name],'char');
fwrite(mesh_file,32*ones(80-6-name_length,1,'int8'),'int8');

% element number sort place
fwrite(mesh_file,0,'uint32');

% write each marker to file
for marker_index=1:length(marker_name_list)
    marker_name=marker_name_list{marker_index};
    if strcmp(marker_name,'point_list') || strcmp(marker_name,'name')
        continue;
    end

    if ~strcmp(grid.(marker_name).type,'stl')
        error('writeMeshSTL: element_type of mesh is not stl format');
    end

    element_list=grid.(marker_name).element_list;

    % write element
    total_element=0;
    element_number=size(element_list,1)/3;
    for element_index=1:element_number
        element=element_list(3*element_index-2:3*element_index,:);

        % write first small element
        point1=element(1,:);
        point2=element(2,:);
        point3=element(3,:);

        d12=point2-point1;
        d23=point3-point2;

        if norm(d12) < eps || norm(d23) < eps
            % small element degeneration to line, discard it
        else
            normal_vector=cross(d12,d23);
            fwrite(mesh_file,normal_vector,'float32');
            fwrite(mesh_file,point1,'float32');
            fwrite(mesh_file,point2,'float32');
            fwrite(mesh_file,point3,'float32');
            fwrite(mesh_file,[0,0],'int8');
            total_element=total_element+1;
        end
            
    end
end

% write element number
fseek(mesh_file,80,-1);
fwrite(mesh_file,total_element,'uint32');

fclose(mesh_file);
clear('mesh_file');

end