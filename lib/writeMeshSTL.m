function writeMeshSTL(filename_mesh,part)
% write Binary STL mesh file
%
% input:
% filename_mesh, part[part.name, part.mesh_list[mesh.element_list, mesh.element_type, mesh.element_ID]]
%
% notice:
% mesh.element_ID if not is 5(tri), will be convert into tri
% part(part.name, part.mesh_list{mesh.element_list, mesh.element_type, mesh.element_ID})
%
if ~isfield(part,'mesh_list')
    error('writeMeshSTL: part lack mesh_list');
end

% check file name
if length(filename_mesh) > 4
    if ~strcmpi(filename_mesh((end-3):end),'.stl')
        filename_mesh=[filename_mesh,'.stl'];
    end
else
    filename_mesh=[filename_mesh,'.stl'];
end

file_mesh=fopen(filename_mesh,'w');

% load part name
if isfield(part,'name')
    part_name=part.name;
else
    part_name=filename_mesh(1:end-4);
end

% write name to file
name_length=length(part_name);
fwrite(file_mesh,['solid ',part_name],'char');
fwrite(file_mesh,32*ones(80-6-name_length,1,'int8'),'int8');

% element number sort place
fwrite(file_mesh,0,'uint32');

% write each element to file
mesh_list=part.mesh_list;
mesh_number=length(mesh_list);
total_element=0;
for mesh_index=1:mesh_number
    mesh=mesh_list{mesh_index};

    % check if exist element_list
    if ~isfield(mesh,'element_list')
        error('writeMeshSTL: part.mesh_list lack element_list');
    end

    element_list=mesh.element_list;

    % write element
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
            fwrite(file_mesh,normal_vector,'float32');
            fwrite(file_mesh,point1,'float32');
            fwrite(file_mesh,point2,'float32');
            fwrite(file_mesh,point3,'float32');
            fwrite(file_mesh,[0,0],'int8');
            total_element=total_element+1;
        end
            
    end
end

% write element number
fseek(file_mesh,80,-1);
fwrite(file_mesh,total_element,'uint32');

fclose(file_mesh);
clear('file_mesh');
end