function [point_list,element_list,marker_list,geometry]=readMeshDataSU2...
    (filename_mesh,scale)
% read mash data from data file
% input mesh_filename(support .su2 file), scale(geometry zoom scale)
% return point_list, element_list(all evaluated marker element), marker_list
%
% point_list is coordinate of all node
% element_list contain element which will be aero function evaluated
% element contain HATSelement
% marker_list contain maker{marker_name,marker_element_number,marker_element}
% marker_element contain contain element(element_type, node_index1, node_index2, ...)
%
if nargin < 2
    scale=[];
end

INFORMATION_FLAG=1;

% cheak file definition
if length(filename_mesh) > 4
    if ~strcmpi(filename_mesh((end-3):end),'.su2')
        filename_mesh=[filename_mesh,'.su2'];
    end
else
    filename_mesh=[filename_mesh,'.su2'];
end
if exist(filename_mesh,'file')==2
    file_mesh=fopen(filename_mesh,'r');
else
    error('readMeshDataSU2: file do not exist')
end

point_list=[];
element_list=[];
marker_list=[];

if isempty(scale)
    scale=1;
end

if INFORMATION_FLAG
    disp('readMeshDataSU2: read mash data begin');
end

dimension_string=fgetl(file_mesh);
dimension_string=strsplit(dimension_string,' ');
dimension=str2double(dimension_string{2});

element_number_string=fgetl(file_mesh);
element_number_string=strsplit(element_number_string,' ');
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
    error('readMeshDataSU2:NPOIN do not exist')
end
point_position=ftell(file_mesh);

% jump over node data first, after read marker, read node data
for point_index=1:point_number
    fgets(file_mesh);
end

% read marker data
marker_number_string=fgetl(file_mesh);
marker_number_string=strsplit(marker_number_string);
if strcmp(marker_number_string{1},'NMARK=')
    marker_number=str2double(marker_number_string{2});
else
    error('readMeshData:NMARK do not exist')
end

marker_list=cell(marker_number,3);
for marker_index=1:marker_number
    % read marker name
    marker_name_string=strsplit(fgetl(file_mesh));
    marker_name=marker_name_string{2};

    % read marker element number
    marker_element_number_string=strsplit(fgetl(file_mesh));
    marker_element_number=str2double(marker_element_number_string{2});

    % read marker element node index
    marker_element=repmat(HATSElement([],[]),marker_element_number,1);
    for element_index=1:marker_element_number
        % new element
        element=HATSElement([],[]);

        marker_element_string=strsplit(fgetl(file_mesh));

        marker_element_type=int8(str2double(marker_element_string{1}));
        element.element_type=marker_element_type;
        switch marker_element_type
            case 3
                element_node_number=2;
            case 5
                element_node_number=3;
            case 9
                element_node_number=4;
        end

        % su2 file point index is start from 0
        element.point_index_list=...
            int32(str2double(marker_element_string(2:1+element_node_number)))+1;
        element.nearby_index_list=int32(zeros(1,element_node_number));

        % give element
        marker_element(element_index)=element;
    end

    marker_list{marker_index,1}=marker_name;
    marker_list{marker_index,2}=marker_element_number;
    marker_list{marker_index,3}=marker_element;
end

% read all point index of marker element
marker_point_number=0;
marker_element_number=0;
for marker_index=1:length(marker_list)
    marker_element_number=marker_element_number+marker_list{marker_index,2};
end
marker_point_index_list=zeros(1,4*marker_element_number);
for marker_index=1:length(marker_list)
    marker_element=marker_list{marker_index,3};
    for element_index=1:marker_list{marker_index,2}
        switch marker_element(element_index).element_type
            case 3
                marker_point_index_list(marker_point_number+1:marker_point_number+2)=...
                    marker_element(element_index).point_index_list;
                marker_point_number=marker_point_number+2;
            case 5
                marker_point_index_list(marker_point_number+1:marker_point_number+3)=...
                    marker_element(element_index).point_index_list;
                marker_point_number=marker_point_number+3;
            case 9
                marker_point_index_list(marker_point_number+1:marker_point_number+4)=...
                    marker_element(element_index).point_index_list;
                marker_point_number=marker_point_number+4;
        end
    end
end

% deleta useless point, get all point index of marker
% mapping_list: old marker_point_index_list to new marker_point_index_list
% read_list: old marker_point_index_list to new marker_point_index_list
% read point should accord with marker_point_index_list
% marker_point_index_list is sort from small to large
marker_point_index_list(marker_point_number+1:end)=[];
[marker_point_index_list,~,mapping_list]=unique(marker_point_index_list);

% updata element_list point index to new list
marker_point_number=0;
for marker_index=1:length(marker_list)
    marker_element=marker_list{marker_index,3};
    for element_index=1:marker_list{marker_index,2}
        switch marker_element(element_index).element_type
            case 2
                element_node_number=2;
                marker_point_number=marker_point_number+2;
            case 5
                element_node_number=3;
                marker_point_number=marker_point_number+3;
            case 9
                element_node_number=4;
                marker_point_number=marker_point_number+4;
        end
        for point_index=1:element_node_number
            marker_element(element_index).point_index_list(point_index)=...
                mapping_list(marker_point_number-element_node_number+point_index);
        end
    end
    marker_list{marker_index,3}=marker_element;
end

marker_point_number=length(marker_point_index_list);

% reread point data
fseek(file_mesh,point_position,"bof");
point_list=zeros(marker_point_number,dimension);
point_index_position=1;
for node_index=1:point_number
    point_string=fgetl(file_mesh);

    % only on marker point will be read
    if marker_point_index_list(point_index_position) == node_index
        point_string=strsplit(point_string,{char(32),char(9)});
        point_list(point_index_position,1:dimension)=str2double(point_string(1:dimension))*scale;
        point_index_position=point_index_position+1;
    end

    if point_index_position > marker_point_number
        % stop reading
        break;
    end
end

geometry.min_bou=min(point_list);
geometry.max_bou=max(point_list);
geometry.dimension=3;

fclose(file_mesh);
clear('file_mesh')

if INFORMATION_FLAG
    disp('readMeshDataSU2: read mash data done!')
end

    function Marker_Element=getMarkerElement(marker_name)
        % return specified marker
        % marker is include all element unit of marker
        % element include element type and point index
        %
        Marker_Element=[];
        for marker_index=1:size(marker_list,1)
            if strcmp(marker_list{marker_index,1},marker_name)
                Marker_Element=marker_list{marker_index,3};
            end
        end
        if isempty(Marker_Element)
            error('getMarkerElement: no marker found')
        end
    end
end