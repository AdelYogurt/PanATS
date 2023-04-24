function user_model=preModelCFG(cfg_filename)
% read model definition from cfg file
%
% copyright Adel 2023.03
%
if nargin < 1
    error('readModelCFG: need input cfg file');
end

if length(cfg_filename) > 4
    if ~strcmpi(cfg_filename((end-3):end),'.cfg')
        cfg_filename=[cfg_filename,'.cfg'];
    end
else
    cfg_filename=[cfg_filename,'.cfg'];
end

% check input file
if exist(cfg_filename,'file') ~= 2
    error('readModelCFG: cfg file do not exist')
end

% initialize model
user_model=struct('g_Point',[],'g_Element',[],'g_Marker',[],...
    'output_inviscid',[],'output_streamline',[],'output_heat',[],'output_viscid',[],'output_FEM',[],...
    'output_post',[]);

cfg_file=fopen(cfg_filename,'r');

% read cfg define line by line
while (~feof(cfg_file))
    string_read=regexprep(fgetl(cfg_file),'\s',''); % read char list and deblank
    if ~isempty(string_read) && string_read(1) ~= '%'
        string_list=strsplit(string_read,{' ','=','{','(',')','}',';',''''});
        parameter=string_list{1};
        value=string_list(2:end);
        if isempty(value{end})
            value=value(1:end-1);
        end

        if isempty(value)
            error('readModelCFG: definition lack value');
        end

        % add parameter, if is number, convert char to number
        digital_value=zeros(1,length(value));
        digital_flag=1;
        for value_index=1:length(value)
            digital_value(value_index)=str2double(value{value_index});
            if isnan(digital_value(value_index))
                digital_flag=0;
                break;
            end
        end

        % if is string
        if digital_flag
            user_model.(parameter)=digital_value;
        else % if is number
            if length(value) == 1
                user_model.(parameter)=value{1};
            else
                user_model.(parameter)=value;
            end
        end
    end
end

fclose(cfg_file);
clear('cfg_file');

% default value
if ~isfield(user_model,'SYMMETRY')
    user_model.SYMMETRY=[];
end
if ~isfield(user_model,'MESH_SCALE')
    user_model.MESH_SCALE=1;
end
if ~isfield(user_model,'MESH_ENCODE')
    user_model.MESH_ENCODE=[];
end
if ~isfield(user_model,'AOA')
    user_model.AOA=0;
end
if ~isfield(user_model,'SIDESLIP_ANGLE')
    user_model.SIDESLIP_ANGLE=0;
end
if ~isfield(user_model,'INFORMATION')
    user_model.INFORMATION=1;
end

% load mesh data
switch user_model.MESH_FORMAT
    case 'SU2'
        % read SU2 format mesh data into mesh format
        [user_model.part_list,user_model.point_list,user_model.geometry]=readMeshSU2...
            (user_model.MESH_FILENAME,user_model.MESH_SCALE,user_model.INFORMATION);

    case 'STL'
        % read each STL format mesh data into STL format
        filename_mesh_list=user_model.MESH_FILENAME;
        if ischar(filename_mesh_list)
            filename_mesh_list={filename_mesh_list};
        end
        part_list=cell(length(filename_mesh_list),1);
        marker_moniter=cell(length(filename_mesh_list),1);

        for file_index=1:length(filename_mesh_list)
            filename_mesh=filename_mesh_list{file_index};
            part_list{file_index}=readMeshSTL...
                (filename_mesh,user_model.MESH_SCALE,user_model.MESH_ENCODE,user_model.INFORMATION);
            marker_moniter{file_index}=part_list{file_index}.name;
        end
        
        % convert STL format into mesh format
        [user_model.part_list,user_model.point_list,user_model.geometry]=convertSTLToMesh(part_list);

        if ~iscell(user_model.part_list)
            user_model.part_list={user_model.part_list};
        end

        % for STL file, if MARKER_MONITORING than will analysis all file
        if ~isfield(user_model,'MARKER_MONITORING')
            user_model.MARKER_MONITORING=marker_moniter;
        end

    case 'INP'
        % read SU2 format mesh data into mesh format
        [user_model.part_list,user_model.point_list,user_model.geometry]=readMeshINP...
            (user_model.MESH_FILENAME,user_model.MESH_SCALE,user_model.INFORMATION);

    case 'LaWGS'
        % read LaWGS format mesh data into mesh format
        part_list=readMeshWGS...
            (user_model.MESH_FILENAME,user_model.MESH_SCALE,user_model.INFORMATION);

        % convert LaWGS format into mesh format
        [user_model.part_list,user_model.point_list,user_model.geometry]=convertWGSToMesh(part_list);
        
end

% generate marker_list
[user_model.marker_list,user_model.element_empty]=preModelMarker(user_model.part_list);

user_model.point_number=size(user_model.point_list,1);
user_model.element_number=sum([user_model.marker_list.element_number]);

if ~isfield(user_model,'MARKER_MONITORING')
    error('preModelCFG: lack marker moniter');
else
    if ischar(user_model.MARKER_MONITORING)
        user_model.MARKER_MONITORING={user_model.MARKER_MONITORING};
    end

    % check MARKER_MONITORING if exist in marker list
    for marker_moniter_index=1:length(user_model.MARKER_MONITORING)
        marker_moniter=user_model.MARKER_MONITORING{marker_moniter_index};
        exist_flag=0;
        for marker_index=1:length(user_model.marker_list)
            if strcmp(user_model.marker_list(marker_index).name,marker_moniter)
                exist_flag=1;
                break;
            end
        end

        if ~exist_flag
            error('preModelCFG: marker monitered do not exist in mesh');
        end
    end
end

end
