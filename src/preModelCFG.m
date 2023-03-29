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
    'inviscid_output',[],'heat_output',[],'viscid_output',[],'FEM_output',[],...
    'post_output',[]);

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
        [user_model.point_list,user_model.element_list,user_model.marker_list,...
            user_model.element_empty,output]=readMeshDataSU2(user_model.MESH_FILENAME,user_model.MESH_SCALE,user_model.INFORMATION);
        user_model.dimension=output.dimension;
        user_model.min_bou=output.min_bou;
        user_model.max_bou=output.max_bou;
    case 'STL'
        [user_model.point_list,user_model.element_list,user_model.marker_list,...
            user_model.element_empty,output,marker_moniter]=readMeshDataSTL(user_model.MESH_FILENAME,user_model.MESH_SCALE,user_model.MESH_ENCODE,user_model.INFORMATION);
        user_model.dimension=output.dimension;
        user_model.min_bou=output.min_bou;
        user_model.max_bou=output.max_bou;

        % for STL file, if MARKER_MONITORING than will analysis all file
        if ~isfield(user_model,'MARKER_MONITORING')
            user_model.MARKER_MONITORING=marker_moniter;
        end
    case 'INP'
        [user_model.point_list,user_model.element_list,user_model.marker_list,...
            user_model.element_empty,output]=readMeshDataINP(user_model.MESH_FILENAME,user_model.MESH_SCALE,user_model.INFORMATION);
        user_model.dimension=output.dimension;
        user_model.min_bou=output.min_bou;
        user_model.max_bou=output.max_bou;
end

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