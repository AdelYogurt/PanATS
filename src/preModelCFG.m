function user_model=preModelCFG(config_filestr)
% read model definition from cfg file
%
% copyright Adel 2023.03
%
if nargin < 1
    error('preModelCFG: need input cfg file');
end

% check input file
if exist(config_filestr,'file') ~= 2
    error('preModelCFG: cfg file do not exist')
end

% initialize model
user_model=struct('element_list',[],'marker_list',[],...
    'output_inviscid',[],'output_streamline',[],'output_boulay',[],'output_viscid',[],'output_heat',[],'output_FEM',[],...
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

end
