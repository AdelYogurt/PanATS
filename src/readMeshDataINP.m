function [point_list,element_list,marker_list,element_empty,geometry]=readMeshDataINP...
    (filename_mesh,scale,INFORMATION)
% read mash data from data file
% input mesh_filename(support .inp file), scale(geometry zoom scale)
% return point_list, element_list, marker_list
%
% point_list is coordinate of all node
% element_list contain element which will be aero function evaluated
% element contain HATSelement
% marker_list is struct list, {marker.name, marker.element_number,marker.element_list}
% marker_element contain contain HATSElement
%
if nargin < 2
    scale=[];
end

% cheak file definition
if length(filename_mesh) > 4
    if ~strcmpi(filename_mesh((end-3):end),'.inp')
        filename_mesh=[filename_mesh,'.inp'];
    end
else
    filename_mesh=[filename_mesh,'.inp'];
end
if exist(filename_mesh,'file')~=2
    error('readMeshDataINP: mesh file do not exist')
else
    file_mesh=fopen(filename_mesh,'r');
end

point_list=[];
element_list=[];
marker_list=[];
element_empty=HATSElement([],[]);

if isempty(scale)
    scale=1;
end

if INFORMATION
    disp('readMeshDataINP: read mash data begin');
end

read_part_flag=0;
read_point_flag=0;
read_element_flag=0;
% loop each line
marker_index=0;
while ~feof(file_mesh)
    string_data=fgetl(file_mesh);
    string_list=strsplit(string_data,{' ',',',char(9)});
    if isempty(string_list{1})
        string_list=string_list(2:end);
    end

    % detact read part
    if strcmp(string_list{1},'*Part')
        read_part_flag=1;

        marker_name=string_list{2}(6:end);
        marker_index=marker_index+1;
        marker_element_list=[];
        continue;
    end

    % detact read point
    if strcmp(string_list{1},'*Node') && read_part_flag
        read_point_flag=1;
        continue;
    end
    % end read point and read element
    if strcmp(string_list{1},'*Element') && read_part_flag
        read_point_flag=0;
        read_element_flag=1;
        continue;
    end

    % end read part
    if strcmp(string_list{1},'*End')
        if ((length(string_list) ==2) && (strcmp(string_list{2},'Part')))
            marker_element_number=size(marker_element_list,1);
            read_part_flag=0;
            read_element_flag=0;
            marker.name=marker_name;
            marker.element_number=marker_element_number;
            marker.element_list=marker_element_list;
            marker_list=[marker_list;marker];
            continue;
        end
    end

    % begin read point data
    if read_point_flag
        point=str2double(string_list(2:4))*scale;
        point_list=[point_list;point];
    end

    % begin read marker element data
    if read_element_flag
        element_node_number=length(string_list)-1;
        switch element_node_number
            case 2
                element_type=int8(3);
            case 3
                element_type=int8(5);
            case 4
                element_type=int8(9);
            otherwise
                error('readMeshDataINP: unknown element type')
        end
        % new element
        element=HATSElement(element_type,int32(str2double(string_list(2:(element_node_number+1)))));
        element.marker_index=marker_index;

        % add element type
        marker_element_list=[marker_element_list;element];
    end

end

fclose(file_mesh);
clear('file_mesh');

geometry.min_bou=min(point_list);
geometry.max_bou=max(point_list);
geometry.dimension=3;

element_list=marker_element_list;

if INFORMATION
    disp('readMeshDataINP: read mash data done');
end

end