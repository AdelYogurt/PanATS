function writeMeshSU2(mesh_filestr,mesh_data,marker_name_list)
% write mesh data into su2 file
%
% input:
% mesh_filestr, mesh_data, marker_name_list(default all markers)
%
% notice:
% mesh_data(single zone): mesh_data.geometry, mesh_data.(marker)
% marker: marker.type, marker.ID, marker.element_list, marker.number_list
% geometry: point_list, dimension
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

mesh_file=fopen(mesh_filestr,'w');

% write dimension
dimension=mesh_data.geometry.dimension;
fprintf(mesh_file,'NDIME= %d\n',dimension);

% write volume element
% check if have marker name 'fluid' or 'fluent' or 'solid'
% if not, check element type to decide which marker is volume element
marker_index=1;
search_flag=1;
marker=[];
while marker_index <= length(marker_name_list)
    marker_name=marker_name_list{marker_index};
    if strcmpi(marker_name,'fluid') || strcmpi(marker_name,'fluent') || strcmpi(marker_name,'solid')
        marker=mesh_data.(marker_name);
        marker_name_list(marker_index)=[];
        search_flag=0;
        break;
    end
end
if search_flag
    % search volume element marker
    for marker_index=1:length(marker_name_list)
        marker_name=marker_name_list{marker_index};
        % temp load marker data
        node_id=mesh_data.(marker_name).ID;
        if node_id == 20 % mixed element
            node_id=mesh_data.(marker_name).element_list(1);
        end

        if (node_id == 5 || node_id == 7)
            if dimension == 2
                marker=mesh_data.(marker_name);
                break;
            end
        elseif node_id == 10 || node_id == 17 || node_id == 14 || node_id == 12
            if dimension == 3
                marker=mesh_data.(marker_name);
                break;
            end
        else
            error('writeMeshSU2: SU2 unsupported element type');
        end
    end
end

if ~isempty(marker)
    % write element number
    if marker.ID == 20 % mixed element
        fprintf(mesh_file,'NELEM= %d\n',sum(marker.number_list));
    else
        fprintf(mesh_file,'NELEM= %d\n',size(marker.element_list,1));
    end

    writeElement(mesh_file,marker);
end

% write point_list
point_list=mesh_data.geometry.point_list;
fprintf(mesh_file,'NPOIN= %d\n',size(point_list,1));
for point_index=1:size(point_list,1)
    fprintf(mesh_file,'%e ',point_list(point_index,1:dimension));
    fprintf(mesh_file,'\n');
end

% write marker element
fprintf(mesh_file,'NMARK= %d\n',length(marker_name_list));
for marker_index=1:length(marker_name_list)
    marker_name=marker_name_list{marker_index};
    marker=mesh_data.(marker_name);

    fprintf(mesh_file,'MARKER_TAG= %s\n',marker_name);
    if marker.ID == 20 % mixed element
        fprintf(mesh_file,'MARKER_ELEMS= %d\n',sum(marker.number_list));
    else
        fprintf(mesh_file,'MARKER_ELEMS= %d\n',size(marker.element_list,1));
    end

    writeElement(mesh_file,marker);
end

fclose(mesh_file);
clear('mesh_file');

    function writeElement(mesh_file,marker)
        ID=marker.ID;
        element_list=marker.element_list;
        number_list=marker.number_list;

        if ID == 20 % mixed element
            data_index=1;
            for element_index=1:length(number_list)
                id=element_list(data_index);
                node_number=number_list(element_index);
                switch id
                    case 3
                        id=3;
                    case 5
                        id=5;
                    case 7
                        id=9;
                    case 10
                        id=10;
                    case 17
                        id=12;
                    case 14
                        id=13;
                    case 12
                        id=14;
                    otherwise
                        error('writeMeshSU2: SU2 unsupported element type');
                end
                fprintf(mesh_file,'%d',id);
                fprintf(mesh_file,' %d',element_list((data_index+1):(data_index+node_number))-1);
                fprintf(mesh_file,'\n');
                data_index=data_index+node_number+1;
            end
        else
            switch ID
                case 3
                    id=3;
                case 5
                    id=5;
                case 7
                    id=9;
                case 10
                    id=10;
                case 17
                    id=12;
                case 14
                    id=13;
                case 12
                    id=14;
                otherwise
                    error('writeMeshSU2: SU2 unsupported element type');
            end
            for element_index=1:size(element_list,1)
                fprintf(mesh_file,'%d',id);
                fprintf(mesh_file,' %d',element_list(element_index,:)-1);
                fprintf(mesh_file,'\n');
            end
        end
    end
end