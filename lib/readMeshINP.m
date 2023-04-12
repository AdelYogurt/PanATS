function [part_list,point_list,geometry]=readMeshINP...
    (filename_mesh,scale,INFORMATION)
% read shell mash data from inp file
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

% cheak filename
if length(filename_mesh) > 4
    if ~strcmpi(filename_mesh((end-3):end),'.inp')
        filename_mesh=[filename_mesh,'.inp'];
    end
else
    filename_mesh=[filename_mesh,'.inp'];
end

% check file if exist
if exist(filename_mesh,'file') ~= 2
    error('readMeshINP: mesh file do not exist')
else
    file_mesh=fopen(filename_mesh,'r');
end

point_list=[];
part_list={};
part_number=0;

if isempty(scale)
    scale=1;
end

if INFORMATION
    disp('readMeshINP: read mash data begin');
end

read_part_flag=0;
read_point_flag=0;
read_element_flag=0;

% loop each line
part_index=0;
point_index_offset=0; % different part point index offset
while ~feof(file_mesh)
    string_data=fgetl(file_mesh);
    string_list=strsplit(string_data,{' ',',',char(9)});
    if isempty(string_list{1})
        string_list=string_list(2:end);
    end

    % detact read part
    if strcmp(string_list{1},'*Part')
        read_part_flag=1;

        % read part
        part_name=string_list{2}(6:end);
        part_index=part_index+1;
        part_mesh_list=[];
        continue;
    end

    % detact read point begain
    if strcmp(string_list{1},'*Node') && read_part_flag
        read_point_flag=1;
        continue;
    end
    % end read point and read element
    if strcmp(string_list{1},'*Element') && read_part_flag && ~read_element_flag
        read_point_flag=0;
        read_element_flag=1;
        element_index=0;

        % begin new mesh read
        mesh_number=0;
        mesh_list={};

        string_list=strsplit(string_list{2},'=');
        mesh_element_type=string_list{2};
        switch mesh_element_type
            case 'S3'
                mesh_element_ID=int8(5);
            case {'S4R','S4'}
                mesh_element_ID=int8(9);
            case 'S8'
                mesh_element_ID=int8(12);
            otherwise
                error('readMeshINP: unknown element type')
        end
        mesh_element_list=[];

        continue;
    end

    % end read part
    if strcmp(string_list{1},'*End')
        if ((length(string_list) ==2) && (strcmp(string_list{2},'Part')))
            read_part_flag=0;
            read_element_flag=0;

            % end mesh read, add old mesh to part mesh list
            mesh.element_list=mesh_element_list;
            mesh.element_type=mesh_element_type;
            mesh.element_ID=mesh_element_ID;

            mesh_number=mesh_number+1;
            mesh_list{mesh_number}=mesh;

            % add old part to part_list
            part.name=part_name;
            part.mesh_list=mesh_list;

            part_number=part_number+1;
            part_list{part_number}=part;

            point_index_offset=int32(size(point_list,1));
            continue;
        end
    end

    % read point data
    if read_point_flag
        point=str2double(string_list(2:4))*scale;
        point_list=[point_list;point];
    end

    % read part element data
    if read_element_flag
        if strcmp(string_list{1},'*Element')
            % detech new mesh, add old mesh to part mesh list
            mesh.element_list=mesh_element_list;
            mesh.element_type=mesh_element_type;
            mesh.element_ID=mesh_element_ID;

            mesh_number=mesh_number+1;
            mesh_list{mesh_number}=mesh;

            % begin new mesh read
            string_list=strsplit(string_list{2},'=');
            mesh_element_type=string_list{2};
            switch mesh_element_type
                case 'S3'
                    mesh_element_ID=int8(5);
                case {'S4R','S4'}
                    mesh_element_ID=int8(9);
                case 'S8'
                    mesh_element_ID=int8(12);
                otherwise
                    error('readMeshINP: unknown element type')
            end
            mesh_element_list=[];

            continue;
        else
            % add element point index
            element_index=element_index+1;

            mesh_element_list=[mesh_element_list;
                int32( str2double( string_list(2:end) ) )+point_index_offset];
        end
    end

end

fclose(file_mesh);
clear('file_mesh');

geometry.min_bou=min(point_list);
geometry.max_bou=max(point_list);
geometry.dimension=3;

if INFORMATION
    disp('readMeshINP: read mash data done');
end

end