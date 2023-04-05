function [point_list,element_list,marker_list,element_empty,geometry,marker_moniter]=readMeshSTL...
    (filename_mesh_list,scale,file_type,INFORMATION)
% read mash data from data file
% input mesh_filename(support .stl file), scale(geometry zoom scale),
% file_type(file encode type, can be empty)
% return point_list, element_list, marker_list
%
% point_list is coordinate of all node
% element_list contain element which will be aero function evaluated
% element contain HATSelement
% marker_list is struct list, {marker.name, marker.element_number,marker.element_list}
% marker_element contain contain element(element_type, node_index1, node_index2, ...)
%
if nargin < 4
    INFORMATION=true(1);
    if nargin < 3
        file_type=[];
        if nargin < 2
            scale=[];
        end
    end
end

geometry_torlance=1e-12;

if ischar(filename_mesh_list)
    filename_mesh_list={filename_mesh_list};
end

point_list=[];
element_list=[];
marker_list=repmat(struct('name',[],'element_number',[],'element_list',[]),length(filename_mesh_list),1);
marker_moniter=cell(length(filename_mesh_list),1); % default marker moniter for STL
element_empty=HATSElement([],[]);

for marker_index=1:length(filename_mesh_list)
    filename_mesh=filename_mesh_list{marker_index};
    % cheak file definition
    if length(filename_mesh) > 4
        if ~strcmpi(filename_mesh((end-3):end),'.stl')
            filename_mesh=[filename_mesh,'.stl'];
        end
    else
        filename_mesh=[filename_mesh,'.stl'];
    end
    if exist(filename_mesh,'file')~=2
        error('readMeshDataSTL: mesh file do not exist')
    end

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

    if INFORMATION
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
        marker_point_list=zeros(marker_element_number*3,3);

        % read element
        repeat_list=[];
        for element_index=1:marker_element_number
            % read normal vector
            vector_normal=fread(file_mesh,3,'float32');

            % read point1 coordinate
            point_1=fread(file_mesh,3,'float32')*scale;
            marker_point_list(element_index*3-2,:)=point_1';

            % read point2 coordinate
            point_2=fread(file_mesh,3,'float32')*scale;
            marker_point_list(element_index*3-1,:)=point_2';

            % read point3 coordinate
            point_3=fread(file_mesh,3,'float32')*scale;
            marker_point_list(element_index*3,:)=point_3';

            attribute=fread(file_mesh,1,'int16');

            if norm(point_1-point_2) < geometry_torlance ||...
                    norm(point_2-point_3) < geometry_torlance ||...
                    norm(point_3-point_1) < geometry_torlance
                repeat_list=[repeat_list;element_index];
            end
        end

        marker_point_list([repeat_list*3-2:repeat_list*3],:)=[];
        marker_element_number=marker_element_number-length(repeat_list);
    else
        % read head
        marker_name=fgetl(file_mesh);
        marker_name=strsplit(marker_name);
        marker_name=marker_name{2};

        % initial sort space
        marker_point_list=zeros(99,3);
        point_number=0;

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

            if norm(point_1-point_2) < geometry_torlance ||...
                    norm(point_2-point_3) < geometry_torlance ||...
                    norm(point_3-point_1) < geometry_torlance
            else
                marker_point_list(point_number+1,:)=point_1';
                marker_point_list(point_number+2,:)=point_2';
                marker_point_list(point_number+3,:)=point_3';
                point_number=point_number+3;
            end

            % add more sort space
            if point_number == size(marker_point_list,1)
                marker_point_list=[marker_point_list;zeros(99,3)];
            end
        end
        marker_point_list=marker_point_list(1:point_number,:);

        marker_element_number=size(marker_point_list,1)/3;
    end

    fclose(file_mesh);
    clear('file_mesh');

    [ADT,index_list]=ADTreePoint(marker_point_list,[],[],geometry_torlance);

    marker_element_list=repmat(HATSElement([],[]),marker_element_number,1);

    geometry.min_bou=ADT.min_bou;
    geometry.max_bou=ADT.max_bou;
    geometry.dimension=3;

    % deleta useless point
    [index_list,~,mapping_list]=unique(index_list);
    % updata element_list point_index to new list
    for element_index=1:marker_element_number
        % new element
        element=HATSElement(int8(5),int32(mapping_list(3*(element_index-1)+(1:3))'));
        element.marker_index=marker_index;
        element.element_index=element_index;

        % give element
        marker_element_list(element_index)=element;
    end

    point_list=[point_list;marker_point_list(index_list,:)];
    marker_list(marker_index).name=marker_name;
    marker_list(marker_index).element_number=marker_element_number;
    marker_list(marker_index).element_list=marker_element_list;

    marker_moniter{marker_index}=marker_name;
end

if INFORMATION
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