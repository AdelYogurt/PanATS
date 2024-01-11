function user_model=preModelCFG(config_filestr)
% read model definition from cfg file
%
% copyright Adel 2023.03
%
if nargin < 1
    error('readModelCFG: need input cfg file');
end

% check input file
if exist(config_filestr,'file') ~= 2
    error('readModelCFG: cfg file do not exist')
end

% initialize model
user_model=struct('element_list',[],'marker_list',[],...
    'output_inviscid',[],'output_streamline',[],'output_boulay',[],'output_heat',[],'output_viscid',[],'output_FEM',[],...
    'output_post',[]);

cfg_file=fopen(config_filestr,'r');

% read cfg define line by line
while (~feof(cfg_file))
    string_read=regexprep(fgetl(cfg_file),{'\s','{','}',' ','(',')',';',''''},''); % read char list and deblank
    if ~isempty(string_read) && string_read(1) ~= '%'
        string_list=strsplit(string_read,{',','='});
        parameter=string_list{1};
        if length(string_list) < 2
            value=[];
        else
            value=string_list(2:end);
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
        % only need marker infomation
        mesh_data=readMeshSU2(user_model.MESH_FILENAME,user_model.MESH_SCALE,false(1),true(1));
    case 'STL'
        mesh_filestr_list=user_model.MESH_FILENAME;
        if ischar(mesh_filestr_list)
            mesh_filestr_list={mesh_filestr_list};
        end
        mesh_data=struct();
        marker_moniter=cell(length(mesh_filestr_list),1);

        for file_index=1:length(mesh_filestr_list)
            mesh_filestr=mesh_filestr_list{file_index};
            mesh_data_new=readMeshSTL(mesh_filestr,user_model.MESH_SCALE,user_model.MESH_ENCODE);
            marker_name=fieldnames(mesh_data_new);marker_name=marker_name{1};
            marker_moniter{file_index}=marker_name;
            mesh_data.(marker_name)=mesh_data_new.(marker_name);
        end
        
        % convert STL format into mesh format
        mesh_data=convertSTLToMesh(mesh_data);

        % for STL file, if MARKER_MONITORING than will analysis all file
        if ~isfield(user_model,'MARKER_MONITORING')
            user_model.MARKER_MONITORING=marker_moniter;
        end
    case 'INP'
        % read SU2 format mesh data into mesh format
        mesh_data=readMeshINP(user_model.MESH_FILENAME,user_model.MESH_SCALE);
    case 'WGS'
        mesh_data=readMeshWGS(user_model.MESH_FILENAME,user_model.MESH_SCALE);

        % convert LaWGS format into mesh format
        mesh_data=convertWGSToMesh(mesh_data);
    case 'CGNS'
        mesh_data=readMeshCGNS(user_model.MESH_FILENAME,user_model.MESH_SCALE);
end

% process marker monitoring
if ~isfield(user_model,'MARKER_MONITORING')
    user_model.MARKER_MONITORING=[];
else
    if isnumeric(user_model.MARKER_MONITORING)
        if user_model.MARKER_MONITORING == 0
            user_model.MARKER_MONITORING=[];
        else
            error('preModelCFG: marker moniter can not be number other than 0');
        end
    end
    
    if ischar(user_model.MARKER_MONITORING)
        user_model.MARKER_MONITORING={user_model.MARKER_MONITORING};
    end

    % check MARKER_MONITORING if exist in mesh data
    for marker_moniter_idx=1:length(user_model.MARKER_MONITORING)
        marker_moniter=user_model.MARKER_MONITORING{marker_moniter_idx};
        exist_flag=0;
        mesh_marker=fieldnames(mesh_data);
        for marker_idx=1:length(mesh_marker)
            if strcmp(mesh_marker{marker_idx},marker_moniter) &&...
                    ~strcmp(mesh_marker{marker_idx},'geometry')
                exist_flag=1;
                break;
            end
        end

        if ~exist_flag
            error('preModelCFG: marker monitered do not exist in mesh');
        end
    end
end

% convert mesh_data to vertex_list
[user_model.element_list,user_model.marker_list]=preModelMesh(mesh_data,user_model.MARKER_MONITORING);
user_model.geometry=mesh_data.geometry;

end
