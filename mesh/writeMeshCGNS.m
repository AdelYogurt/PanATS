function writeMeshCGNS(mesh_filestr,mesh_data,marker_name_list)
% write mesh data into cgns file
% need cgns4m package
%
% input:
% mesh_filestr, mesh_data, marker_name_list(default all markers)
%
% notice:
% mesh_data(single zone): mesh_data.name, mesh_data.point_list, mesh_data.(marker)
% marker: marker.type, marker.ID, marker.element_list, marker.number_list
% notice: marker which has the same name of file is volume element
%
if nargin < 3
    marker_name_list=[];
end
if isempty(marker_name_list)
    marker_name_list=fieldnames(mesh_data);
end
marker_index=1;
while marker_index <= length(marker_name_list)
    marker_name=marker_name_list{marker_index};
    if strcmp(marker_name,'geometry')
        marker_name_list(marker_index)=[];
    else
        marker_index=marker_index+1;
    end
end

[~,mesh_filename,~]=fileparts(mesh_filestr);

% open the CGNS file.
ierr=cg_set_file_type(CG_FILE_ADF); chk_error(ierr);
[mesh_file, ierr]=cg_open(mesh_filestr, CG_MODE_WRITE); chk_error(ierr);

% point list
point_list=mesh_data.geometry.point_list;
dimension=size(point_list, 2);

if isfield(mesh_data,mesh_filename) || isfield(mesh_data,upper(mesh_filename))
    if ~isfield(mesh_data,mesh_filename)
        mesh_filename=upper(mesh_filename);
    end
    marker=mesh_data.(mesh_filename);
    element_dimension=getDimension(size(marker.element_list,2),marker.type);
    if strcmp(marker.type,'MIXED2') || strcmp(marker.type,'MIXED3')
        element_number=length(marker.number_list);
    else
        element_number=size(marker.element_list,1);
    end
else
    element_dimension=dimension;
    element_number=0;
%     % load all marker to find out max dimension
%     dimension_list=zeros(length(marker_name_list),1);
%     for marker_index=1:length(marker_name_list)
%         marker=mesh_data.(marker_name_list{marker_index});
%         dimension_list(marker_index)=getDimension(size(marker.element_list,2),marker.type);
%     end
%     [element_dimension,marker_index]=max(dimension_list);
%     marker=mesh_data.(marker_name_list{marker_index});
%     if strcmp(marker.type,'MIXED2') || strcmp(marker.type,'MIXED3')
%         element_number=length(marker.number_list);
%     else
%         element_number=size(marker.element_list,1);
%     end
end

% Create base
[base_index, ierr]=cg_base_write(mesh_file, 'Base', element_dimension, dimension); chk_error(ierr);

% Number of vertices and elements
data_size=[size(point_list, 1), element_number, zeros(1, 7)];
% Create zone
[zone_index, ierr]=cg_zone_write(mesh_file, base_index, mesh_filename, data_size, ...
    Unstructured); chk_error(ierr);

% Write mesh_data coordinates (must use SIDS-standard names here)
[~, ierr]=cg_coord_write(mesh_file, base_index, zone_index, ...
    RealDouble, 'CoordinateX', point_list(:, 1)); chk_error(ierr);
[~, ierr]=cg_coord_write(mesh_file, base_index, zone_index, ...
    RealDouble, 'CoordinateY', point_list(:, 2)); chk_error(ierr);
if dimension == 3
    [~, ierr]=cg_coord_write(mesh_file, base_index, zone_index, ...
        RealDouble, 'CoordinateZ', point_list(:, 3)); chk_error(ierr);
end

offset=0;
for marker_index=1:length(marker_name_list)
    marker_name=marker_name_list{marker_index};
    marker=mesh_data.(marker_name);

    ID=marker.ID;
    element_list=marker.element_list;

    if ID == 20
        element_number=length(marker.number_list);
    else
        element_number=size(element_list,1);
    end

    % Write element connectivity. We must permute elems, but we don't need to
    % cast the data type to integer explicitly (MEX function does it for us).
    [~, ierr]=cg_section_write(mesh_file, base_index, zone_index, marker_name, ...
        ID, offset+1, offset+element_number, 0, element_list'); chk_error(ierr);

    offset=offset+size(element_list,1);

end

ierr=cg_close(mesh_file); chk_error(ierr);
clear('mesh_file');

end

function element_dimension=getDimension(node_number, typestr)
% Obtain the element-type ID and dimension of elements
switch (node_number)
    case 1
        if strcmpi(typestr, 'MIXED2')
            element_dimension=2;
        elseif strcmpi(typestr, 'MIXED3')
            element_dimension=3;
        end
    case 2
        element_dimension=1;
    case 3
        if ~isempty(typestr) && upper(typestr(1)) == 'B'
            element_dimension=1;
        else
            element_dimension=2;
        end
    case 4
        if ~isempty(typestr) && upper(typestr(1)) == 'HF'
            element_dimension=2;
        else
            element_dimension=3;
        end
    case 5
        element_dimension=3;
    case 6
        if (~isempty(typestr) && upper(typestr(1)) == 'P')
            element_dimension=3;
        else
            element_dimension=2;
        end
    case 8
        if ~isempty(typestr) && upper(typestr(1)) == 'HF'
            element_dimension=2;
        else
            element_dimension=3;
        end
    case 9
        element_dimension=2;
    case {10,13,14,15,18,20,27}
        element_dimension=3;
    otherwise
        error('ERROR: unknown element type with %d nodes.', node_number);
end
end

function chk_error(ierr)
% Check whether CGNS returned an error code. If so, get error message
if ierr
    error(['Error: ', cg_get_error()]);
end
end
