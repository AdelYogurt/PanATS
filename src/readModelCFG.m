function my_model=readModelCFG(cfg_filename)
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
my_model=struct('g_Point',[],'g_Element',[],'g_Marker',[],...
    'geometry_output',[],'inviscid_output',[],'heat_output',[],'FEM_output',[],...
    'post',[]);

cfg_file=fopen(cfg_filename,'r');

% read cfg define line by line
while (~feof(cfg_file))
    string_read=fgetl(cfg_file);
    if ~isempty(string_read) && string_read(1) ~= '%'
        string_list=strsplit(string_read,{' ','=','{','(',')','}',';'});
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
            my_model.(parameter)=value;
        else % if is number
            my_model.(parameter)=digital_value;
        end
    end
end

fclose(cfg_file);
clear('cfg_file');
end