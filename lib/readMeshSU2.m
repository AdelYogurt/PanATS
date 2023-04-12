function [part_list,point_list,geometry]=readMeshSU2...
    (filename_mesh,scale,INFORMATION)
% read mash data from su2 file, only include marker mesh data
%
% input:
% filename_mesh(support .inp file), scale(geometry zoom scale), ...
% INFORMATION(true or false)
%
% output:
% point_list, part_list
%
% notice:
% point_list is coordinate of all node
% part_list{part.name, part.mesh_list{mesh.element_list, mesh.element_type, mesh.element_ID}}
%
if nargin < 3
    INFORMATION=true(1);
    if nargin < 2
        scale=[];
    end
end


% cheak file definition
if length(filename_mesh) > 4
    if ~strcmpi(filename_mesh((end-3):end),'.su2')
        filename_mesh=[filename_mesh,'.su2'];
    end
else
    filename_mesh=[filename_mesh,'.su2'];
end

% check file if exist
if exist(filename_mesh,'file') ~= 2
    error('readMeshSU2: mesh file do not exist')
else
    file_mesh=fopen(filename_mesh,'r');
end

point_list=[];
part_list={};

if isempty(scale)
    scale=1;
end

if INFORMATION
    disp('readMeshSU2: read mash data begin');
end

% read mesh dimension
dimension_string=strsplit(fgetl(file_mesh),' ');
dimension=str2double(dimension_string{2});

% read volume element number
element_number_string=strsplit(fgetl(file_mesh),' ');
element_number=str2double(element_number_string{2});

% jump over volume element data
for element_index=1:element_number
    fgets(file_mesh);
end

% read node data
point_number_string=fgetl(file_mesh);
point_number_string=strsplit(point_number_string,' ');
if strcmp(point_number_string{1},'NPOIN=')
    point_number=str2double(point_number_string{2});
else
    error('readMeshSU2:NPOIN do not exist')
end
point_position=ftell(file_mesh);

% jump over node data first, after read marker, read node data
for point_index=1:point_number
    fgets(file_mesh);
end

% read marker data
part_number_string=strsplit(fgetl(file_mesh));
if strcmp(part_number_string{1},'NMARK=')
    part_number=str2double(part_number_string{2});
else
    error('readMeshData:NMARK do not exist')
end

part_point_number=0;
part_list=cell(part_number,1);
for part_index=1:part_number
    % read marker name
    part_name_string=strsplit(fgetl(file_mesh));
    part_name=part_name_string{2};

    % read marker element number
    part_element_number_string=strsplit(fgetl(file_mesh));
    part_element_number=str2double(part_element_number_string{2});

    % initialize mesh list
    mesh_tri.element_ID=5;
    mesh_tri.element_type='S3';
    mesh_tri.element_list=zeros(part_element_number,3);
    mesh_tri.element_number=0;

    mesh_quad.element_ID=9;
    mesh_quad.element_ID='S4R';
    mesh_quad.element_list=zeros(part_element_number,4);
    mesh_quad.element_number=0;

    mesh_list={mesh_tri,mesh_quad};

    % read marker element node index
    for element_index=1:part_element_number
        element_string=strsplit(fgetl(file_mesh));
        element_ID=int8(str2double(element_string{1}));
        
        switch element_ID
            case 5
                element_node_number=3;
                mesh_index=1;
            case 9
                element_node_number=4;
                mesh_index=2;
        end

        % su2 file point index is start from 0, so we need to add 1
        point_index=int32(str2double(element_string(2:1+element_node_number)))+1;

        % give element
        addMesh(mesh_index,point_index);
    end

    % delete empty mesh
    mesh_index=1;
    while mesh_index <= length(mesh_list)
        if mesh_list{mesh_index}.element_number == 0
            mesh_list(mesh_index)=[];
        else
            element_number=mesh_list{mesh_index}.element_number;
            mesh_list{mesh_index}.element_list=...
                mesh_list{mesh_index}.element_list(1:element_number,:);
            part_point_number=part_point_number+element_number*element_node_number;
            mesh_index=mesh_index+1;
        end
    end

    part.name=part_name;
    part.mesh_list=mesh_list;

    part_list{part_index}=part;
end

% read all point index of part element
flag_number=0;
total_point_index_list=zeros(1,part_point_number);
for part_index=1:length(part_list)
    mesh_list=part_list{part_index}.mesh_list;

    for mesh_index=1:length(mesh_list)
        mesh=mesh_list{mesh_index};

        element_number=size(mesh.element_list,1);

        switch mesh.element_ID
            case 5
                element_node_number=3;
            case 9
                element_node_number=4;
        end
        
        mesh_point_number=element_number*element_node_number;

        % add all point_index to total_point_index_list
        total_point_index_list((flag_number+1):(flag_number+mesh_point_number))=...
            reshape(mesh.element_list',1,mesh_point_number);
        flag_number=flag_number+mesh_point_number;
    end
end

% deleta useless point, get all point index of marker
% mapping_list: old total_point_index_list to new total_point_index_list
% total_point_index_list is sorted from small to large
[total_point_index_list,~,mapping_list]=unique(total_point_index_list);

% updata element_list point index to new list
flag_number=0;
for part_index=1:length(part_list)
    mesh_list=part_list{part_index}.mesh_list;

    for mesh_index=1:length(mesh_list)
        mesh=mesh_list{mesh_index};

        element_number=size(mesh.element_list,1);

        switch mesh.element_ID
            case 5
                element_node_number=3;
            case 9
                element_node_number=4;
        end
        
        mesh_point_number=element_number*element_node_number;

        % generate new element list
        mesh.element_list=reshape(mapping_list((flag_number+1):(flag_number+mesh_point_number))...
            ,element_node_number,element_number)';

        flag_number=flag_number+mesh_point_number;

        mesh_list{mesh_index}=mesh;
    end

    part_list{part_index}.mesh_list=mesh_list;
end

% reread point data
part_point_number=length(total_point_index_list);
fseek(file_mesh,point_position,"bof");
point_list=zeros(part_point_number,3);
flag_point=1;
for node_index=1:point_number
    point_string=fgetl(file_mesh);

    % only on marker point will be read
    if total_point_index_list(flag_point) == node_index

        point_string=strsplit(point_string,{char(32),char(9)});
        point_list(flag_point,1:3)=str2double(point_string(1:3))*scale;
        flag_point=flag_point+1;
    end

    if flag_point > part_point_number
        % stop reading
        break;
    end
end

fclose(file_mesh);
clear('file_mesh')

geometry.min_bou=min(point_list);
geometry.max_bou=max(point_list);
geometry.dimension=3;

if INFORMATION
    disp('readMeshSU2: read mash data done!')
end

    function addMesh(mesh_index,point_index)
        mesh_list{mesh_index}.element_number=mesh_list{mesh_index}.element_number+1;
        mesh_list{mesh_index}.element_list(mesh_list{mesh_index}.element_number,:)=point_index;
    end

end
