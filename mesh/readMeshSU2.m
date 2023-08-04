function [grid,geometry]=readMeshSU2(mesh_filestr,scale,READ_VOLUME,ONLY_MARKER)
% read mesh data from su2 file
%
% input:
% filename_mesh(support .su2 file), scale(geometry zoom scale), ...
% READ_VOLUME(true or false), ONLY_MARKER(only marker data, true or false)
%
% output:
% grid, geometry
%
% notice:
% INDEX OF POINT HAVE PLUS ONE!
% grid(single zone): grid.name, grid.point_list, grid.(marker)
% marker: marker.type, marker.element_list, marker.number_list, marker.ID_list
% notice: marker which has the same name of file is volume element
%
if nargin < 4
    ONLY_MARKER=[];
    if nargin < 3
        READ_VOLUME=[];
        if nargin < 2
            scale=[];
        end
    end
end

if isempty(ONLY_MARKER), ONLY_MARKER=false(1);end
if isempty(READ_VOLUME), READ_VOLUME=true(1);end

% cheak file definition
if numel(mesh_filestr) > 4
    if ~strcmpi(mesh_filestr((end-3):end),'.su2')
        mesh_filestr=[mesh_filestr,'.su2'];
    end
else
    mesh_filestr=[mesh_filestr,'.su2'];
end

% check file if exist
if exist(mesh_filestr,'file') ~= 2
    error('readMeshSU2: mesh file do not exist')
else
    mesh_file=fopen(mesh_filestr,'r');
end

[~,mesh_filename,~]=fileparts(mesh_filestr);

if isempty(scale)
    scale=1;
end

while ~feof(mesh_file)

    % string=strsplit(fgetl(mesh_file),{char(9),char(32)});
    string=strsplit(fgetl(mesh_file));

    if strcmp(string{1},'NDIME=')
        dimension=str2double(string{2});
    elseif strcmp(string{1},'NPOIN=')
        % read node data
        point_number=str2double(string{2});

        if ONLY_MARKER
            point_position=ftell(mesh_file);
            % jump over node data first, after read marker, read node data
            for point_index=1:point_number
                fgets(mesh_file);
            end
        else
            point_list=zeros(point_number,dimension);
            for point_index=1:point_number
                node_string=strsplit(fgetl(mesh_file));
                for dimension_index=1:dimension
                    point_list(point_index,dimension_index)=str2double(node_string{dimension_index})*scale;
                end
            end
            grid.point_list=point_list;
        end
    elseif strcmp(string{1},'NELEM=')
        % read volume element number
        element_number=str2double(string{2});
        if READ_VOLUME && ~ONLY_MARKER
            % read volume element
            [type,element_list,number_list]=readElement(element_number,mesh_file);
            grid.(mesh_filename).type=type;
            grid.(mesh_filename).element_list=element_list;
            grid.(mesh_filename).number_list=number_list;
        else
            % jump over volume element data
            for index=1:element_number
                fgetl(mesh_file);
            end
        end
    elseif strcmp(string{1},'NMARK=')
        % read marker data
        marker_number=str2double(string{2});
        marker_name_list=cell(marker_number,1);
        for marker_index=1:marker_number
            % read marker name
            marker_name_string=strsplit(fgetl(mesh_file));
            marker_name=regexprep(marker_name_string{2},{'[',']','"','''','-'},'');
            marker_name_list{marker_index}=marker_name;

            % read marker element number
            marker_element_number_string=strsplit(fgetl(mesh_file));
            marker_element_number=str2double(marker_element_number_string{2});

            % read marker element
            [type,element_list,number_list]=readElement(marker_element_number,mesh_file);

            grid.(marker_name).type=type;
            grid.(marker_name).element_list=element_list;
            grid.(marker_name).number_list=number_list;
        end
    end

end

if ONLY_MARKER
    % read all point index of part element
    total_point_index_list=[];
    for marker_index=1:numel(marker_name_list)
        marker_name=marker_name_list{marker_index};
        total_point_index_list=[total_point_index_list;grid.(marker_name).element_list(:)];
    end

    % deleta useless point, get all point index of marker
    % map_list: old total_point_index_list to new total_point_index_list
    % total_point_index_list is sorted from small to large
    [total_point_index_list,~,map_list]=unique(total_point_index_list);

    % updata element_list point index to new list
    flag_number=0;
    for marker_index=1:numel(marker_name_list)
        marker_name=marker_name_list{marker_index};
        marker_node_number=numel(grid.(marker_name).element_list);

        % generate new element list
        grid.(marker_name).element_list(1:marker_node_number)=map_list((flag_number+1):(flag_number+marker_node_number));
        
        flag_number=flag_number+marker_node_number;
    end

    % reread point data
    point_number=numel(total_point_index_list);
    fseek(mesh_file,point_position,"bof");
    point_list=zeros(point_number,dimension);
    flag_point=1;
    for point_index=1:point_number
        point_string=fgetl(mesh_file);

        % only on marker point will be read
        if total_point_index_list(flag_point) == point_index
            % point_string=strsplit(point_string,{char(32),char(9)});
            point_string=strsplit(point_string);
            for dimension_index=1:dimension
                point_list(flag_point,dimension_index)=str2double(point_string{dimension_index})*scale;
            end
            flag_point=flag_point+1;
        end

        if flag_point > point_number
            % stop reading
            break;
        end
    end
    grid.point_list=point_list;
end

fclose(mesh_file);
clear('mesh_file')

grid.name=mesh_filename;
geometry.dimension=dimension;

    function [type,element_list,number_list,ID_list]=readElement(element_number,mesh_file)
        % read element_number element from mesh_file
        %
        SAME_TYPE=true(1);

        add_number=element_number*2;
        element_list=zeros(add_number,1,'uint32');
        number_list=zeros(element_number,1,'uint8');
        ID_list=zeros(element_number,1,'uint8');

        % read marker element node index
        node_index=0;
        for element_index=1:element_number
            element_string=strsplit(fgetl(mesh_file));
            identifier=str2double(element_string{1});
            node_number=idSU2Number(identifier);
            number_list(element_index,1)=idSU2Number(identifier);
            ID_list(element_index,1)=idSU2ID(identifier);

            % append element
            if node_index+node_number > numel(element_list)
                element_list=[element_list;zeros(add_number,1,'uint32');];
            end

            for node_idx=1:node_number
                element_list(node_index+node_idx)=str2double(element_string{1+node_idx})+1;
            end
            node_index=node_index+node_number;
        end
        element_list=element_list(1:node_index);

        % check if have same type
        node_number=number_list(1);
        for element_index=1:element_number
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

    end

end

function node_number=idSU2Number(identifier)
switch identifier
    case 3
        node_number=2;
    case 5
        node_number=3;
    case 9
        node_number=4;
    case 10
        node_number=4;
    case 12
        node_number=8;
    case 13
        node_number=6;
    case 14
        node_number=5;
    otherwise
        error('idSU2Number: unknown identifier')
end

end

function id=idSU2ID(identifier)
switch identifier
    case 3
        id=3;
    case 5
        id=5;
    case 9
        id=7;
    case 10
        id=10;
    case 12
        id=17;
    case 13
        id=14;
    case 14
        id=12;
    otherwise
        error('idSU2Number: unknown identifier')
end

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
