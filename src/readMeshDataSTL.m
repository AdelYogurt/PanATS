function [point_list,element_list,marker_list,geometry]=readMeshDataSTL...
    (filename_mesh,scale,file_type)
% read mash data from data file
% input mesh_filename(support .stl file), scale(geometry zoom scale),
% file_type(file encode type, can be empty)
% return point_list, element_list, marker_list
%
% point_list is coordinate of all node
% element_list contain element which will be aero function evaluated
% element contain element_type, node_index1, node_index2, ...
% marker_list contain maker{marker_name,marker_element_number,marker_element}
% marker_element contain contain element(element_type, node_index1, node_index2, ...)
%
if nargin < 3
    file_type=[];
    if nargin < 2
        scale=[];
    end
end

INFORMATION_FLAG=1;

% cheak file definition
if length(filename_mesh) > 4
    if ~strcmpi(filename_mesh((end-3):end),'.stl')
        filename_mesh=[filename_mesh,'.stl'];
    end
else
    filename_mesh=[filename_mesh,'.stl'];
end
if exist(filename_mesh,'file')~=2
    error('readMeshDataSTL: file do not exist')
end

point_list=[];
element_list=[];
marker_list=[];

torlance=1e-12;

if isempty(scale)
    scale=1;
end

% detech file type
if isempty(file_type)
    file_mesh=fopen(filename_mesh,'rb');
    string_head=fread(file_mesh,80,'int8');
    
    if string_head(end)==32
        file_type='binary';
        fseek(file_mesh,0,'bof');
    else
        file_type='ASCII';
        fclose(file_mesh);
        clear('file_mesh');
        
        % reopen as ASCII file
        file_mesh=fopen(filename_mesh,'r');
    end
end

if INFORMATION_FLAG
    disp(['readMeshDataSTL: read mash data begin, file type read as ',file_type]);
end

if strcmp(file_type,'binary')
    % read head
    marker_name=fread(file_mesh,80,'int8');
    marker_name=char(marker_name');
    marker_name=strsplit(marker_name);
    marker_name=marker_name{2};
    
    % read total element_number
    marker_element_number=fread(file_mesh,1,'int32');
    point_list=zeros(marker_element_number*3,3);
    
    % read element
    repeat_list=[];
    for element_index=1:marker_element_number
        % read normal vector
        vector_normal=fread(file_mesh,3,'float32');
        
        % read point1 coordinate
        point_1=fread(file_mesh,3,'float32')*scale;
        point_list(element_index*3-2,:)=point_1';
        
        % read point2 coordinate
        point_2=fread(file_mesh,3,'float32')*scale;
        point_list(element_index*3-1,:)=point_2';
        
        % read point3 coordinate
        point_3=fread(file_mesh,3,'float32')*scale;
        point_list(element_index*3,:)=point_3';
        
        attribute=fread(file_mesh,1,'int16');
        
        if norm(point_1-point_2) < torlance ||...
                norm(point_2-point_3) < torlance ||...
                norm(point_3-point_1) < torlance
            repeat_list=[repeat_list;element_index];
        end
    end
    
    point_list([repeat_list*3-2:repeat_list*3],:)=[];
    marker_element_number=marker_element_number-length(repeat_list);
else
    % read head
    marker_name=fgetl(file_mesh);
    marker_name=strsplit(marker_name);
    marker_name=marker_name{2};
    
    % read normal vector
    vector_normal_string=strsplit(fgetl(file_mesh));
    while ~strcmp(vector_normal_string{1},'endsolid')
        fgetl(file_mesh);
        
        % read point1 coordinate
        point_1=getASCIIPoint(file_mesh)*scale;
        
        % read point2 coordinate
        point_2=getASCIIPoint(file_mesh)*scale;
        
        % read point3 coordinate
        point_3=getASCIIPoint(file_mesh)*scale;
        
        fgetl(file_mesh);fgetl(file_mesh);
        
        % read normal vector
        vector_normal_string=strsplit(fgetl(file_mesh));
        
        if norm(point_1-point_2) < torlance ||...
                norm(point_2-point_3) < torlance ||...
                norm(point_3-point_1) < torlance
        else
            point_list=[point_list;point_1'];
            point_list=[point_list;point_2'];
            point_list=[point_list;point_3'];
        end
    end
    
    marker_element_number=size(point_list,1)/3;
end

fclose(file_mesh);
clear('file_mesh');

[ADT,index_list]=ADTreePoint(point_list,[],[],torlance);
marker_element=zeros(marker_element_number,4);

geometry.min_bou=ADT.min_bou;
geometry.max_bou=ADT.max_bou;
geometry.dimension=3;

% deleta useless point
[index_list,~,mapping_list]=unique(index_list);
% updata element_list point index to new list
for element_index=1:marker_element_number
    marker_element(element_index,1)=5;
    for point_index=1:3
        marker_element(element_index,1+point_index)=...
            mapping_list(3*(element_index-1)+point_index);
    end
end

point_list=point_list(index_list,:);
element_list=marker_element;
marker_list={marker_name,marker_element_number,marker_element};

if INFORMATION_FLAG
    disp('readMeshDataSTL: read mash data done!');
end

    function point=getASCIIPoint(file_mesh)
        % read ASCII point data
        %
        point_string=strsplit(fgetl(file_mesh));
        point=zeros(3,1);
        for node_index=1:3
            point(node_index)=str2double(point_string{2+node_index});
        end
    end
end