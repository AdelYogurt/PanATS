function mesh_data=readMeshSU2(mesh_filestr,scale,READ_VOLUME,ONLY_MARKER)
% read mesh data from su2 file
%
% input:
% mesh_filestr(support .su2 file), scale(geometry zoom scale), ...
% READ_VOLUME(true or false), ONLY_MARKER(only marker data, true or false)
%
% output:
% mesh_data, geometry
%
% notice:
% INDEX OF POINT HAVE PLUS ONE!
% mesh_data(single zone): mesh_data.geometry, mesh_data.(marker)
% marker: marker.type, marker.ID, marker.element_list, marker.number_list
% geometry: point_list, dimension
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

        point_list=textscan(mesh_file,'%f');
        point_list=reshape(point_list{1},[],point_number)';
        point_list=point_list(:,1:dimension);

        if scale ~= 1
            point_list=point_list*scale;
        end
    elseif strcmp(string{1},'NELEM=')
        % read volume element number
        element_number=str2double(string{2});
        if READ_VOLUME && ~ONLY_MARKER
            % read volume element
            [type,ID,element_list,number_list]=readElement(element_number,mesh_file);
            mesh_data.(mesh_filename).type=type;
            mesh_data.(mesh_filename).ID=ID;
            mesh_data.(mesh_filename).element_list=element_list;
            mesh_data.(mesh_filename).number_list=number_list;
        else
            textscan(mesh_file,'%d');
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
            [type,ID,element_list,number_list]=readElement(marker_element_number,mesh_file);

            mesh_data.(marker_name).type=type;
            mesh_data.(marker_name).ID=ID;
            mesh_data.(marker_name).element_list=element_list;
            mesh_data.(marker_name).number_list=number_list;
        end
    end

end

if ONLY_MARKER
    % read all point index of part element
    total_point_index_list=[];
    for marker_index=1:length(marker_name_list)
        marker_name=marker_name_list{marker_index};
        total_point_index_list=[total_point_index_list;mesh_data.(marker_name).element_list(:)];
    end

    % deleta useless point, get all point index of marker
    % map_list: old total_point_index_list to new total_point_index_list
    % total_point_index_list is sorted from small to large
    [total_point_index_list,~,map_list]=unique(total_point_index_list);

    % updata element_list point index to new list
    flag_number=0;
    for marker_index=1:length(marker_name_list)
        marker_name=marker_name_list{marker_index};
        marker_node_number=numel(mesh_data.(marker_name).element_list);

        % generate new element list
        mesh_data.(marker_name).element_list(1:marker_node_number)=map_list((flag_number+1):(flag_number+marker_node_number));

        flag_number=flag_number+marker_node_number;
    end

    % delete no marker point data
    point_list=point_list(total_point_index_list,:);
end

fclose(mesh_file);
clear('mesh_file')

geometry.dimension=dimension;
geometry.point_list=point_list;
mesh_data.geometry=geometry;

    function [type,ID,element_list,number_list]=readElement(element_number,mesh_file)
        % read element_number element from mesh_file
        %
        SAME_TYPE=true(1);

        % read unstructure element list
        element_list=textscan(mesh_file,'%d');
        element_list=element_list{1}+1;

        % judge if have element index in each line
        switch element_list(1)-1
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
                error('idSU2Number: unknown SU2 identifier')
        end
        if element_list(node_number+2) == 1
            offset=1;
        else
            offset=0;
        end

        % change SU2 ID into ID and create number_list
        number_list=zeros(element_number,1);
        index_list=zeros(element_number,1);
        data_number=length(element_list);
        data_index=1;
        element_index=1;
        while data_index < data_number
            switch element_list(data_index)-1
                case 3
                    node_number=2;
                    id=3;
                case 5
                    node_number=3;
                    id=5;
                case 9
                    node_number=4;
                    id=7;
                case 10
                    node_number=4;
                    id=10;
                case 12
                    node_number=8;
                    id=17;
                case 13
                    node_number=6;
                    id=14;
                case 14
                    node_number=5;
                    id=12;
                otherwise
                    error('idSU2Number: unknown SU2 identifier')
            end
            element_list(data_index)=id;
            number_list(element_index)=node_number;

            data_index=data_index+1+node_number;
            if offset
                index_list(element_index)=data_index;
                data_index=data_index+1;
            end
            element_index=element_index+1;
        end

        % remove element index
        if offset
            element_list(index_list)=[];
        end

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
            ID=element_list(1);
            element_list=reshape(element_list,node_number+1,[])';
            element_list=element_list(:,2:node_number+1);
            number_list=number_list(1);
        else
            ID=20;
        end
        type=convertIDToType(ID,dimension);
    end
end

function type=convertIDToType(id,dimension)
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
    case 20
        if dimension == 2
            type = 'MIXED2';
        elseif dimension == 3
            type = 'MIXED3';
        else
            error('convertIDToType: for mixed meshes, dimension must be 2 or 3.');
        end
    otherwise
        error('convertIDToType: unknown identifier')
end

end
