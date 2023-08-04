function [grid,geometry]=readMeshINP(mesh_filestr,scale)
% read mesh data from inp file
%
% input:
% filename_mesh(support .inp file), scale(geometry zoom scale)
%
% output:
% grid, geometry
%
% notice:
% point_list is coordinate of all node
% grid(single zone): grid.name, grid.point_list, grid.(marker)
% marker: marker.type, marker.element_list, marker.number_list, marker.ID_list
% notice: marker which has the same name of file is volume element
%
if nargin < 3
    INFORMATION=true(1);
    if nargin < 2
        scale=[];
    end
end

% cheak filename
if length(mesh_filestr) > 4
    if ~strcmpi(mesh_filestr((end-3):end),'.inp')
        mesh_filestr=[mesh_filestr,'.inp'];
    end
else
    mesh_filestr=[mesh_filestr,'.inp'];
end

% check file if exist
if exist(mesh_filestr,'file') ~= 2
    error('readMeshINP: mesh file do not exist')
else
    mesh_file=fopen(mesh_filestr,'r');
end

[~,mesh_filename,~]=fileparts(mesh_filestr);

point_list=[];
dimension=3;
add_number=50;

if isempty(scale)
    scale=1;
end

read_marker_flag=0;
read_point_flag=0;
read_element_flag=0;

% loop each line
point_index_offset=0; % different part point index offset
while ~feof(mesh_file)
    string_data=fgetl(mesh_file);
    string_list=strsplit(string_data,{' ',',',char(9)});
    if isempty(string_list{1})
        string_list=string_list(2:end);
    end

    if strcmp(string_list{1},'*Part')
        % detact read part
        read_marker_flag=1;

        % read part
        marker_name=regexprep(string_list{2}(6:end),{'[',']','"','''','-'},'');

        node_index=0;
        element_index=0;

        element_list=zeros(add_number,1);
        number_list=zeros(add_number,1);
        ID_list=zeros(add_number,1);
    elseif strcmp(string_list{1},'*Node') && read_marker_flag
        % detact read point begain
        read_point_flag=1;

    elseif strcmp(string_list{1},'*Element') && read_marker_flag
        % end read point and read element
        read_point_flag=0;
        read_element_flag=1;

        string_list=strsplit(string_list{2},'=');
        mesh_element_type=string_list{2};
        switch mesh_element_type
            case 'S3'
                element_id=uint8(5);
                node_number=3;
            case {'S4R','S4'}
                element_id=uint8(7);
                node_number=4;
            case 'S8'
                element_id=uint8(17);
                node_number=8;
            otherwise
                error('readMeshINP: unknown element type')
        end

    elseif strcmp(string_list{1},'*End')
        % end read part
        if ((length(string_list) ==2) && (strcmp(string_list{2},'Part')))
            read_marker_flag=0;
            read_element_flag=0;

            element_list=element_list(1:node_index);
            number_list=number_list(1:element_index);
            ID_list=ID_list(1:element_index);

            SAME_TYPE=true(1);

            % check if have same type
            node_number=number_list(1);
            for element_index=1:numel(number_list)
                if number_list(element_index) ~= node_number
                    SAME_TYPE=false(1);
                    break;
                end
            end

             % if same type, reshape element
             if SAME_TYPE
                 element_list=reshape(element_list,node_number,[])';
                 number_list=number_list(1);
                 ID_list=ID_list(1);
                 type=idType(ID_list);
             else
                 if dimension == 2
                     type = 'MIXED2';
                 elseif dimension == 3
                     type = 'MIXED3';
                 else
                     error('readMeshSU2.readElement: for mixed meshes, dimension must be 2 or 3.');
                 end
             end

            % end mesh read, add marker
            grid.(marker_name).type=type;
            grid.(marker_name).element_list=element_list;
            grid.(marker_name).number_list=number_list;
            grid.(marker_name).ID_list=ID_list;

            point_index_offset=uint32(size(point_list,1));
        end
    else
        % read point data
        if read_point_flag
            point=str2double(string_list(2:1+dimension))*scale;
            point_list=[point_list;point];
        end

        % read part element data
        if read_element_flag
            % add element point index
            if numel(element_list) < node_index+node_number
                element_list=[element_list;zeros(add_number,1)];
            end

            for index=1:node_number
                element_list(node_index+index)=uint32(str2double(string_list(1+index)))+point_index_offset;
            end

            if numel(ID_list) < element_index+1
                ID_list=[ID_list;zeros(add_number,1)];
                number_list=[number_list;zeros(add_number,1)];
            end

            ID_list(element_index+1)=element_id;
            number_list(element_index+1)=node_number;

            node_index=node_index+node_number;
            element_index=element_index+1;
        end
    end
end

fclose(mesh_file);
clear('mesh_file');

grid.point_list=point_list;

grid.name=mesh_filename;
geometry.dimension=dimension;

end

function type=idType(id)
switch id
    case 3
        type='BAR_2';
    case 5
        type='TRI_3';
    case 7
        type='QUAD_4';
    case 10
        type='TETRA_4';
    case 17
        type='HEXA_8';
    case 14
        type='PENTA_6';
    case 12
        type='PYRA_5';
    otherwise
        error('idType: unknown identifier')
end

end

