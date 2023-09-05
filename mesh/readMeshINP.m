function mesh_data=readMeshINP(mesh_filestr,scale)
% read mesh data from inp file
%
% input:
% mesh_filestr(support .inp file), scale(geometry zoom scale)
%
% output:
% mesh_data, geometry
%
% notice:
% point_list is coordinate of all node
% mesh_data(single zone): mesh_data.geometry, mesh_data.(marker)
% marker: marker.type, marker.ID, marker.element_list, marker.number_list
% geometry: point_list, dimension
%
if nargin < 2
    scale=[];
end

% check file if exist
if exist(mesh_filestr,'file') ~= 2
    error('readMeshINP: mesh file do not exist')
else
    mesh_file=fopen(mesh_filestr,'r');
end

[~,mesh_filename,~]=fileparts(mesh_filestr);

point_list=[];
dimension=3;

if isempty(scale)
    scale=1;
end

% loop each line
point_index_offset=0; % different part point index offset
while ~feof(mesh_file)
    string_data=fgetl(mesh_file);
    string_list=strsplit(string_data,{' ',',',char(9)});

    if strcmp(string_list{1},'*Part')
        % read part
        marker_name=regexprep(string_list{2}(6:end),{'[',']','"','''','-'},'');

        [point_list_new,type,ID,element_list,number_list,...
            Nset,Elset,material]=readPart(point_index_offset,mesh_file);

        % end mesh read, add marker
        mesh_data.(marker_name).type=type;
        mesh_data.(marker_name).ID=ID;
        mesh_data.(marker_name).element_list=element_list;
        mesh_data.(marker_name).number_list=number_list;
        mesh_data.(marker_name).Nset=Nset;
        mesh_data.(marker_name).Elset=Elset;
        mesh_data.(marker_name).material=material;

        point_list=[point_list;point_list_new];
        point_index_offset=int32(size(point_list,1));
    end
    
    if strcmp(string_list{1},'*Assembly')
        break;
    end
end

% residual_data=textscan(mesh_file,'%s');
% residual_data=[string(string_data);residual_data{1}];

residual_data=textscan(mesh_file,'%s','delimiter','\n','whitespace','');
residual_data=[string(string_data);residual_data{1}];

fclose(mesh_file);
clear('mesh_file');

if scale ~= 1
    point_list=point_list*scale;
end

geometry.dimension=dimension;
geometry.point_list=point_list;
geometry.residual_data=residual_data;
mesh_data.geometry=geometry;

    function [point_list,type,ID,element_list,number_list,...
            Nset,Elset,material]=readPart(point_index_offset,mesh_file)
        % read part data
        %

        % read point list
        str_data=fgetl(mesh_file);
        str_list=strsplit(str_data,{' ',',',char(9)});
        if strcmp(str_list{1},'*Node')
            point_list=textscan(mesh_file,'%f,%f,%f,%f');
            point_list=[point_list{2:4}];
        end

        str_data=fgetl(mesh_file);
        str_list=strsplit(str_data,{' ',',',char(9)});
        element_list=[];
        Nset=[];
        Elset=[];
        material=[];
        while ~strcmp(str_list{1},'*End')

            % read element list
            if strcmp(str_list{1},'*Element')
                % element type
                str_list=strsplit(str_list{2},'=');
                type_new=str_list{2};
                [ID_new,number_list_new]=convertTypeToID(type_new);

                % read element list
                scan_format=['%d',repmat(',%d',1,number_list_new)];
                element_list_new=textscan(mesh_file,scan_format);
                element_list_new=[element_list_new{2:(1+number_list_new)}]+(point_index_offset);

                if isempty(element_list)
                    element_list=element_list_new;
                    type=type_new;
                    ID=ID_new;
                    number_list=number_list_new;
                else
                    % mixed element
                    number_list_new=repmat(number_list_new,size(element_list_new,1),1);

                    if number_list_new(end) ~= number_list(end) && length(number_list) == 1
                        number_list=repmat(number_list,size(element_list,1),1);
                        element_list=[repmat(ID,size(element_list,1),1),element_list]';
                        element_list=element_list(:);
                        type='MIXED3';
                    end

                    element_list_new=[repmat(ID_new,size(element_list_new,1),1),element_list_new]';
                    element_list=[element_list;element_list_new(:)];
                    ID=20;
                    number_list=[number_list;number_list_new];
                end
            end

            % read node setting
            if strcmp(str_list{1},'*Nset')
                % name
                temp_list=strsplit(str_list{2},'=');
                set_name=regexprep(temp_list{2},{'[',']','"','-'},'');
                node_index=textscan(mesh_file,'%d,');
                node_index=node_index{1};
                if length(str_list) > 2 && strcmp(str_list{end},'generate')
                    % mean generate from start to end
                    node_index=(node_index(1):node_index(3):node_index(2))';
                end
                Nset.(set_name)=node_index;
            end

            % read element setting
            if strcmp(str_list{1},'*Elset')
                % name
                temp_list=strsplit(str_list{2},'=');
                set_name=regexprep(temp_list{2},{'[',']','"','-'},'');
                elememt_index=textscan(mesh_file,'%d,');
                elememt_index=elememt_index{1};
                if length(str_list) > 2 && strcmp(str_list{end},'generate')
                    % mean generate from start to end
                    elememt_index=(elememt_index(1):elememt_index(3):elememt_index(2))';
                end
                Elset.(set_name)=elememt_index;
            end

            % read materials setting
            if strcmp(str_list{1},'*Solid')
                temp_list=strsplit(str_list{3},'=');
                Nset_name=regexprep(temp_list{2},{'[',']','"','-'},'');
                temp_list=strsplit(str_list{4},'=');
                marterial_name=regexprep(temp_list{2},{'[',']','"','-'},'');
                material.(Nset_name)=marterial_name;
            end

            str_data=fgetl(mesh_file);
            str_list=strsplit(str_data,{' ',',',char(9)});
        end

    end
end

function [ID,node_number]=convertTypeToID(type)
% inp version of ID and type converter
%
switch type
    case 'S3'
        ID=uint8(5);
        node_number=3;
    case {'S4R','S4'}
        ID=uint8(7);
        node_number=4;
    case {'S8R','S8R5'}
        ID=uint8(8);
        node_number=8;
    case 'S9R'
        ID=uint8(9);
        node_number=9;
    case 'C3D4'
        ID=uint8(10);
        node_number=4;
    case {'C3D8','C3D8R'}
        ID=uint8(17);
        node_number=8;
    otherwise
        error('readMeshINP: unknown element type')
end
end
