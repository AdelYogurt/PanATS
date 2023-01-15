function [geometry,point_list,element_list,marker_list]=readMeshDataSU2...
    (filename_mesh,evaluate_marker_name_list)
% read mash data from data file
% input mesh_filename, support .su2 file
% return point_list, element_list, marker_list

% point_list is coordinate of all node
% element_list contain element(element_type, node_index1, node_index2, ...)
% marker_list contain make{marker_name,marker_element_number,marker_element} 
% marker_element contain contain element(element_type, node_index1, node_index2, ...)
%
disp('readMeshData: read mash data begin')

point_list=[];
element_list=[];
marker_list=[];

if ~strcmp(filename_mesh((end-3):end),'.su2')
    filename_mesh=[filename_mesh,'.su2'];
end

if exist(filename_mesh,'file')==2
    file_mesh=fopen(filename_mesh,'r');
else
    error('readMeshData: file do not exist')
end

dimension_string=fgetl(file_mesh);
dimension_string=strsplit(dimension_string);
dimension=str2double(dimension_string{2});

element_number_string=fgetl(file_mesh);
element_number_string=strsplit(element_number_string);
element_number=str2double(element_number_string{2});

% jump over volume element data
for rank_index=1:element_number
    fgetl(file_mesh);
end
 
% read node data
node_number_string=fgetl(file_mesh);
node_number_string=strsplit(node_number_string);
if strcmp(node_number_string{1},'NPOIN=')
    node_number=str2double(node_number_string{2});
else
    error('readMeshData:NPOIN do not exist')
end
point_list=zeros(node_number,dimension);
for node_index=1:node_number
    node_string=strsplit(fgetl(file_mesh));
    if strcmp(node_string{1},'')
       node_string=node_string(2:end); 
    end
    for dimension_index=1:dimension
        point_list(node_index,dimension_index)=str2double(node_string{dimension_index});
    end
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
    if dimension == 2
        marker_element=zeros(marker_element_number,3);
    elseif dimension == 3
        marker_element=zeros(marker_element_number,5);
    end
    for marker_element_index=1:marker_element_number
        marker_element_string=strsplit(fgetl(file_mesh));
        marker_element_type=str2double(marker_element_string{1});
        marker_element(marker_element_index,1)=marker_element_type;
        switch marker_element_type
            case 3
                marker_element_index_number=2;
            case 5
                marker_element_index_number=3;
            case 9
                marker_element_index_number=4;
        end
        for node_index=1:marker_element_index_number
            marker_element(marker_element_index,node_index+1)=str2double(marker_element_string{node_index+1});
        end
    end
    
    marker_list{marker_index,1}=marker_name;
    marker_list{marker_index,2}=marker_element_number;
    marker_list{marker_index,3}=marker_element;
end
fclose(file_mesh);
clear('file_mesh')

geometry.min_bou=min(point_list);
geometry.max_bou=max(point_list);
geometry.dimension=3;

% get all element of marker and create HATSElemet list g_Marker_Element
element_list=[]; % all element needed to be evaluate list, pointer to HATSElemet
for marker_name_index=1:size(evaluate_marker_name_list,2)
    element_list=[element_list;getMarkerElement(evaluate_marker_name_list{marker_name_index})];
end

% read all point on marker
marker_point_number=0;
marker_point_index_list=zeros(1,4*size(element_list,1));
for element_index=1:size(element_list,1)
    element_type=element_list(element_index,1);
    switch element_type
        case 5
            marker_point_index_list(marker_point_number+1:marker_point_number+3)=...
                element_list(element_index,2:4);
            marker_point_number=marker_point_number+3;
        case 9
            marker_point_index_list(marker_point_number+1:marker_point_number+4)=...
                element_list(element_index,2:5);
            marker_point_number=marker_point_number+4;
    end
end

% deleta useless point
marker_point_index_list(marker_point_number+1:end)=[];
[marker_point_index_list,~,mapping_list]=unique(marker_point_index_list);

% updata element_list point index to new list
marker_point_number=0;
for element_index=1:size(element_list,1)
    element_type=element_list(element_index,1);
    switch element_type
        case 5
            point_number=3;
            marker_point_number=marker_point_number+3;
        case 9
            point_number=4;
            marker_point_number=marker_point_number+4;
    end
    for point_index=1:point_number
        element_list(element_index,1+point_index)=...
            mapping_list(marker_point_number-point_number+point_index);
    end
end

point_list=point_list(marker_point_index_list+1,:);

disp('readMeshData: read mash data done!')

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