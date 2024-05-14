function mesh_data=readMeshSTL(mesh_filestr,scale,file_type,geometry_torlance)
% read mesh data from stl file
%
% input:
% mesh_filestr(support .stl file), scale(geometry zoom scale), ...
% file_type(stl code format, binary/ASCII), ...
% geometry_torlance(default is 1e-12)
%
% output:
% mesh_data
%
% notice:
% mesh_data(single zone): mesh_data.(marker)
% marker: marker.type, marker.element_list
%
if nargin < 4
    geometry_torlance=1e-12;
    if nargin < 3
        file_type=[];
        if nargin < 2
            scale=[];
        end
    end
end

if isempty(scale)
    scale=1;
end

% check file if exist
if exist(mesh_filestr,'file')~=2
    error('readMeshSTL: mesh file do not exist')
end

dimension=3;

% detech file type
if isempty(file_type)
    mesh_file=fopen(mesh_filestr,'rb');
    string_head=fread(mesh_file,80,'int8');

    if string_head(end)==32 || string_head(end)==0
        file_type='binary';
        fseek(mesh_file,0,'bof');
    else
        file_type='ASCII';
        fclose(mesh_file);
        clear('mesh_file');

        % reopen as ASCII file
        mesh_file=fopen(mesh_filestr,'r');
    end
end

[~,marker_name,~]=fileparts(mesh_filestr);

if strcmp(file_type,'binary')
    % read head
    marker_name_string=fread(mesh_file,80,'int8');

    % read total element_number
    element_number=fread(mesh_file,1,'int32');
    element_list=zeros(3*element_number,3);

    % read element
    overlap_list=[];
    for element_index=1:element_number
        % read normal vector
        vector_normal=fread(mesh_file,3,'float32');

        % read point coordinate
        point_1=fread(mesh_file,3,'float32');
        point_2=fread(mesh_file,3,'float32');
        point_3=fread(mesh_file,3,'float32');

        element_list(3*element_index-2:3*element_index,:)=...
            [point_1,point_2,point_3]';

        attribute=fread(mesh_file,1,'int16');

        if norm(point_1-point_2) < geometry_torlance ||...
                norm(point_2-point_3) < geometry_torlance ||...
                norm(point_3-point_1) < geometry_torlance
            overlap_list=[overlap_list;element_index];
        end
    end

    element_list([3*overlap_list-2;3*overlap_list-1;overlap_list],:)=[];
    element_number=element_number-length(overlap_list);
else
    % read head
    marker_name_string=fgetl(mesh_file);

    % initial sort space
    element_list=zeros(99,3);
    element_number=0;

    % read normal vector
    vector_normal_string=strsplit(fgetl(mesh_file));
    while ~strcmp(vector_normal_string{1},'endsolid')
        fgetl(mesh_file);

        % read point1 coordinate
        point_1=getASCIIPoint(mesh_file);
        point_2=getASCIIPoint(mesh_file);
        point_3=getASCIIPoint(mesh_file);

        fgetl(mesh_file);fgetl(mesh_file);

        % read normal vector
        vector_normal_string=strsplit(fgetl(mesh_file));

        if norm(point_1-point_2) < geometry_torlance ||...
                norm(point_2-point_3) < geometry_torlance ||...
                norm(point_3-point_1) < geometry_torlance
        else
            element_number=element_number+1;
            element_list(3*element_number-2:3*element_number,:)=...
                [point_1,point_2,point_3]';
        end

        % add more sort space
        if element_number*3 == size(element_list,1)
            element_list=[element_list;zeros(99,3)];
        end
    end
    element_list=element_list(1:element_number*3,:);

    element_number=size(element_list,1)/3;
end

fclose(mesh_file);
clear('mesh_file');

if scale ~= 1
    element_list=element_list*scale;
end

mesh_data.(marker_name).element_list=element_list;
mesh_data.(marker_name).type='stl';

    function point=getASCIIPoint(mesh_file)
        % read ASCII point data
        %
        point_string=strsplit(fgetl(mesh_file));
        point=zeros(3,1);
        for node_idx=1:3
            point(node_idx)=str2double(point_string{2+node_idx});
        end
    end
end