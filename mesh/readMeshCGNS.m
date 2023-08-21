function mesh_data=readMeshCGNS(mesh_filestr,scale)
% read mesh data from cgns file
% need cgns4m package
%
% input:
% filename_mesh(support .su2 file), scale(geometry zoom scale)
%
% output:
% mesh_data, geometry
%
% notice:
% mesh_data(single zone): mesh_data.geometry, mesh_data.(marker)
% marker: marker.type, marker.ID, marker.element_list, marker.number_list
% geometry: point_list, dimension
%
if nargin < 2
    scale=[];
end

% check file if exist
if exist(mesh_filestr,'file') ~= 2
    error('readMeshCGNS: mesh file do not exist')
end
[~,mesh_filename,~]=fileparts(mesh_filestr);
if isempty(scale)
    scale=1;
end

% only support one base one zone now
base_index=1;
zone_index=1;

% open file
[mesh_file, ierr]=cg_open(mesh_filestr, CG_MODE_READ); chk_error(ierr);

% read point
[basename, cell_dimension, dimension, ierr]=cg_base_read(mesh_file, base_index); chk_error(ierr);
[zonename, size, ierr]=cg_zone_read(mesh_file, base_index, zone_index); chk_error(ierr);

% Define the range of vertices
rmin=1; % lower range idx of vertices
rmax=size(1); % upper range idx of vertices

% Read mesh_data coordinates (must use SIDS-standard names here)
point_list=zeros(rmax(1), dimension);
[point_list(:, 1), ierr]=cg_coord_read(mesh_file, base_index, zone_index, ...
    'CoordinateX', RealDouble, rmin, rmax, point_list(:, 1)); chk_error(ierr);
[point_list(:, 2), ierr]=cg_coord_read(mesh_file, base_index, zone_index, ...
    'CoordinateY', RealDouble, rmin, rmax, point_list(:, 2)); chk_error(ierr);
if (dimension == 3)
    [point_list(:, 3), ierr]=cg_coord_read(mesh_file, base_index, zone_index, ...
        'CoordinateZ', RealDouble, rmin, rmax, point_list(:, 3)); chk_error(ierr);
end
if scale ~= 1
    point_list=point_list*scale;
end

% read section
[marker_number, ierr]=cg_nsections(mesh_file, base_index, zone_index); chk_error(ierr);
for marker_index=1:marker_number

    [marker_name, ID, istart, iend, nbndry, pflag, ierr]=cg_section_read(mesh_file, ...
        base_index, zone_index, marker_index); chk_error(ierr);

    [node_number, type]=convertIDToType(ID, cell_dimension);

    % Get element connectivity
    [node_number_total, ierr]=cg_ElementDataSize(mesh_file, base_index, zone_index, ...
        marker_index); chk_error(ierr);
    num_elems=node_number_total / node_number; % number of elements (except for MIXED)
    element_list=zeros(node_number, num_elems); % Element connectivity is permuted in CGNS
    parent_data=[];

    [element_list, parent_data, ierr]=cg_elements_read(mesh_file, base_index, zone_index, ...
        marker_index, element_list, parent_data); chk_error(ierr);
    element_list=element_list'; % Permute the connectivity back

    if ID == 20
        % for mixed elemenet, separate into number_list and ID_list
        [element_list,number_list]=getNumberElement(element_list);
    else
        number_list=node_number;
    end

    mesh_data.(marker_name).type=type;
    mesh_data.(marker_name).ID=ID;
    mesh_data.(marker_name).element_list=element_list;
    mesh_data.(marker_name).number_list=number_list;
end

ierr=cg_close(mesh_file); chk_error(ierr);
clear('mesh_file');

geometry.point_list=point_list;
geometry.dimension=dimension;

mesh_data.geometry=geometry;

end

function [node_number, type]=convertIDToType(ID, cell_dimension)
% Obtain a string of element type
switch (ID)
    case NODE
        type='NODE';
        node_number=1;
    case BAR_2
        type='BAR_2';
        node_number=2;
    case BAR_3
        type='BAR_3';
        node_number=3;
    case TRI_3
        type='TRI_3';
        node_number=3;
    case TRI_6
        type='TRI_6';
        node_number=6;
    case QUAD_4
        type='QUAD_4';
        node_number=4;
    case QUAD_8
        type='QUAD_8';
        node_number=8;
    case QUAD_9
        type='QUAD_9';
        node_number=9;
    case TETRA_4
        type='TETRA_4';
        node_number=4;
    case TETRA_10
        type='TETRA_10';
        node_number=10;
    case PYRA_5
        type='PYRA_5';
        node_number=5;
    case PYRA_13
        type='PYRA_13';
        node_number=13;
    case PYRA_14
        type='PYRA_14';
        node_number=14;
    case PENTA_6
        type='PENTA_6';
        node_number=6;
    case PENTA_15
        type='PENTA_15';
        node_number=15;
    case PENTA_18
        type='PENTA_18';
        node_number=18;
    case HEXA_8
        type='HEXA_8';
        node_number=8;
    case HEXA_20
        type='HEXA_20';
        node_number=20;
    case HEXA_27
        type='HEXA_27';
        node_number=27;
    case MIXED
        if cell_dimension == 2
            type='MIXED2';
        elseif cell_dimension == 3
            type='MIXED3';
        else
            error('For mixed meshes, dimension must be 2 or 3.');
        end
        node_number=1;

    otherwise
        error('Error: unknown element type');
end
end

function [element_list,number_list]=getNumberElement(element_list)
% obtain number list by element list
%

% calculate element number
number_list=zeros(length(element_list),1,'int64');
element_index=ones(1,1,'int64');

index=1;
data_number=length(element_list);
while index < data_number
    switch element_list(index)
        case 5
            number_list(element_index)=3;
        case 7
            number_list(element_index)=4;
        case 6
            number_list(element_index)=6;
        case 8
            number_list(element_index)=8;
        case 9
            number_list(element_index)=9;
        case 10
            number_list(element_index)=4;
        case 12
            number_list(element_index)=5;
        case 14
            number_list(element_index)=6;
        case 17
            number_list(element_index)=8;
        case 11
            number_list(element_index)=10;
        case 21
            number_list(element_index)=13;
        case 13
            number_list(element_index)=14;
        case 15
            number_list(element_index)=15;
        case 16
            number_list(element_index)=18;
        case 18
            number_list(element_index)=20;
        case 19
            number_list(element_index)=27;
        otherwise
            error('ERROR: unknown element type in elems(%d).', element_index);
    end

    index=index+number_list(element_index)+1;
    element_index=element_index+1;
end

element_number=element_index-1;
number_list=number_list(1:element_number);
end

function chk_error(ierr)
% Check whether CGNS returned an error code. If so, get error message
if ierr
    error(['Error: ', cg_get_error()]);
end
end
