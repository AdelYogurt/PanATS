function [point_list,element_list,marker_list,geometry]=readMeshDataINP...
    (filename_mesh,scale)
% read mash data from data file
% input mesh_filename(support .inp file), scale(geometry zoom scale)
% return point_list, element_list, marker_list
%
% point_list is coordinate of all node
% element_list contain element which will be aero function evaluated
% element contain element_type, node_index1, node_index2, ...
% marker_list contain maker{marker_name,marker_element_number,marker_element}
% marker_element contain contain element(element_type, node_index1, node_index2, ...)
%
if nargin < 2
    scale=[];
end

INFORMATION_FLAG=1;

% cheak file definition
if length(filename_mesh) > 4
    if ~strcmpi(filename_mesh((end-3):end),'.inp')
        filename_mesh=[filename_mesh,'.inp'];
    end
else
    filename_mesh=[filename_mesh,'.inp'];
end
if exist(filename_mesh,'file')~=2
    error('readMeshDataINP: file do not exist')
else
    file_mesh=fopen(filename_mesh,'r');
end

point_list=[];
element_list=[];
marker_list=[];

if isempty(scale)
    scale=1;
end

if INFORMATION_FLAG
    disp('readMeshDataINP: read mash data begin');
end

read_part_flag=0;
read_point_flag=0;
read_element_flag=0;
% loop each line
while ~feof(file_mesh)
    string_data=fgetl(file_mesh);
    string_list=strsplit(string_data,{' ',',',char(9)});
    if isempty(string_list{1})
        string_list=string_list(2:end);
    end

    % detact read part
    if strcmp(string_list{1},'*Part')
        marker_name=string_list{2}(6:end);
        read_part_flag=1;

        marker_element=[];
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
            marker_element_number=size(marker_element,1);
            read_part_flag=0;
            read_element_flag=0;
            marker={marker_name,marker_element_number,marker_element};
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
        element_point_number=length(string_list)-1;
        element=zeros(1,element_point_number+1);

        % add node index
        element(2:(element_point_number+1))=str2double(string_list(2:(element_point_number+1)));

        switch element_point_number
            case 2
                element_type=3;
            case 3
                element_type=5;
            case 4
                element_type=9;
            otherwise
                error('readMeshDataINP: unknown element type')
        end

        % add element type
        element(1)=element_type;
        marker_element=[marker_element;element];
    end

end

fclose(file_mesh);
clear('file_mesh');

geometry.min_bou=min(point_list);
geometry.max_bou=max(point_list);
geometry.dimension=3;

element_list=marker_element;

if INFORMATION_FLAG
    disp('readMeshDataINP: read mash data done');
end

end