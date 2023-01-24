function g_model=preModelCFG(cfg_filename)
% read model definition from cfg file
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
g_model=struct('g_Point',[],'g_Element',[],'g_Marker',[],...
    'inviscid_output',[],'heat_output',[],'viscid_output',[],'FEM_output',[],...
    'post_output',[]);

cfg_file=fopen(cfg_filename,'r');

% read cfg define line by line
while (~feof(cfg_file))
    string_read=fgetl(cfg_file);
    if ~isempty(string_read) && string_read(1) ~= '%'
        string_list=strsplit(string_read,{' ','=','{','(',')','}',';',''''});
        parameter=string_list{1};
        value=string_list{2:end};

        %         % check out if exist empty
        %         value_index=1;
        %         while (value_index <= length(value))
        %             if isempty(value{value_index})
        %                 value(value_index=[]);
        %             end
        %             value_index=value_index+1;
        %         end
        if isempty(value)
            error('readModelCFG: definition lack value');
        end

        % add parameter
        digital_value=str2double(value);
        if isnan(digital_value(1)) % if is string
            g_model.(parameter)=value;
        else % if is number
            g_model.(parameter)=digital_value;
        end
    end
end

fclose(cfg_file);
clear('cfg_file');

% default value
if ~isfield(g_model,'SYMMETRY')
    g_model.SYMMETRY=[];
end
if ~isfield(g_model,'MESH_SCALE')
    g_model.MESH_SCALE=1;
end
if ~isfield(g_model,'MESH_ENCODE')
    g_model.MESH_ENCODE=[];
end
if ~isfield(g_model,'AOA')
    g_model.AOA=0;
end
if ~isfield(g_model,'SIDESLIP_ANGLE')
    g_model.SIDESLIP_ANGLE=0;
end

% load mesh data
switch g_model.MESH_FORMAT
    case 'SU2'
        [g_model.point_list,g_model.element_list,g_model.marker_list,...
            output]=readMeshDataSU2(g_model.MESH_FILENAME,g_model.MESH_SCALE);
        g_model.dimension=output.dimension;
        g_model.min_bou=output.min_bou;
        g_model.max_bou=output.max_bou;
    case 'STL'
        [g_model.point_list,g_model.element_list,g_model.marker_list,...
            output,marker_moniter]=readMeshDataSTL(g_model.MESH_FILENAME,g_model.MESH_SCALE,g_model.MESH_ENCODE);
        g_model.dimension=output.dimension;
        g_model.min_bou=output.min_bou;
        g_model.max_bou=output.max_bou;

        % for STL file, if MARKER_MONITORING than will analysis all file
        if ~isfield(g_model,'MARKER_MONITORING')
            g_model.MARKER_MONITORING=marker_moniter;
        end
    case 'INP'
        [g_model.point_list,g_model.element_list,g_model.marker_list,...
            output]=readMeshDataINP(g_model.MESH_FILENAME,g_model.MESH_SCALE);
        g_model.dimension=output.dimension;
        g_model.min_bou=output.min_bou;
        g_model.max_bou=output.max_bou;
end

if ~isfield(g_model,'MARKER_MONITORING')
    error('preModelCFG: lack marker moniter');
else
    if ischar(g_model.MARKER_MONITORING)
        g_model.MARKER_MONITORING={g_model.MARKER_MONITORING};
    end
end

end