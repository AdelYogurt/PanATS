classdef PanATSConfig < handle
    % dict of SU2 config
    properties
        config_filedir;
        config_filename;
        dump_filename='config.cfg';
        data_dict=struct();
    end

    % main function
    methods
        function self=PanATSConfig(config_filestr)
            % initialize config, process input config_filestr
            % read config parameter from config file
            %
            if nargin < 1
                error('PanATSConfig: need input cfg file');
            end

            if length(config_filestr) > 4
                if ~strcmpi(config_filestr((end-3):end),'.cfg')
                    config_filestr=[config_filestr,'.cfg'];
                end
            else
                config_filestr=[config_filestr,'.cfg'];
            end

            % check if input file exist
            if ~exist(config_filestr,'file')
                error('PanATSConfig: cfg file do not exist');
            end

            [self.config_filedir,temp_filename,~]=fileparts(config_filestr);
            self.config_filename=[temp_filename,'.cfg'];

            self.setDefault();
            self.configRead();
        end

        function dump(self,cfg_filestr)
            if nargin < 2
                cfg_filestr=self.dump_filename;
            end
            self.configWrite(cfg_filestr);
        end

        function configRead(self,config_filestr)
            if nargin < 2
                config_filestr=fullfile(self.config_filedir,self.config_filename);
            else
                if length(config_filestr) > 4
                    if ~strcmpi(config_filestr((end-3):end),'.cfg')
                        config_filestr=[config_filestr,'.cfg'];
                    end
                else
                    config_filestr=[config_filestr,'.cfg'];
                end

                % check if input file exist
                if ~exist(config_filestr,'file')
                    error('PanATSConfig: cfg file do not exist');
                end

                [self.config_filedir,temp_filename,~]=fileparts(which(config_filestr));
                self.config_filename=[temp_filename,'.cfg'];
            end

            cfg_file=fopen(config_filestr,'r');

            % read cfg define line by line
            while (~feof(cfg_file))
                string_read=regexprep(fgetl(cfg_file),'\s',''); % read char list and deblank
                if ~isempty(string_read) && string_read(1) ~= '%'
                    % if end with \, means need to read next line
                    while string_read(end) == '\'
                        string_read=string_read(1:end-1);
                        string_temp=regexprep(fgetl(cfg_file),'\s',''); % read char list and deblank
                        if len(string_temp.split('=')) > 1
                            error('PanATSConfig: statement found after end')
                        end
                        if string_read(1) ~= '%'
                            string_read=[string_read,string_temp];
                        end
                    end

                    % replace {}();' by ,
                    % and split string by , and =
                    string_read=regexprep(string_read,'[{}();''|]',',');
                    string_list=strsplit(string_read,{',','='});
                    parameter=string_list{1};
                    value=string_list(2:end);
                    if isempty(value{end})
                        value=value(1:end-1);
                    end

                    if isempty(value),value=[];end

                    % add parameter, if is number, convert char to number
                    for value_idx=1:length(value)
                        digital_value=str2double(value{value_idx});
                        if ~isnan(digital_value)
                            value{value_idx}=digital_value;
                        end
                    end

                    if length(value) == 1
                        self.data_dict.(parameter)=value{1};
                    else
                        self.data_dict.(parameter)=value;
                    end
                end
            end

            fclose(cfg_file);
            clear('cfg_file');

        end

        function configWrite(self,cfg_filestr)
            % write data_dict to cfg file
            %
            cfg_file=fopen(cfg_filestr,'w');
            parameter_list=fieldnames(self.data_dict);
            
            % write each parameter
            for param_idx=1:length(parameter_list)
                % printf parameter
                parameter=parameter_list{param_idx};
                fprintf(cfg_file,'%s=',parameter);

                % printf value
                % notice if value more than one, use ()
                value=self.data_dict.(parameter);
                if ~iscell(value), value={value}; end

                if length(value)>1 || contains(parameter,'MARKER')
                    fprintf(cfg_file,'(');
                    for value_idx=1:length(value)-1
                        printValue(value{value_idx})
                        fprintf(cfg_file,',');
                    end
                    printValue(value{end})
                    fprintf(cfg_file,')');
                else
                    printValue(value{1})
                end

                fprintf(cfg_file,'\n');
            end

            fclose(cfg_file);
            clear('cfg_file');

            function printValue(value)
                if isnumeric(value)
                    fprintf(cfg_file,'%s',num2str(value));
                else
                    fprintf(cfg_file,'%s',value);
                end
            end
        end

        function setDefault(self)
            % set default value to config
            %
            self.data_dict.SOLVER='PANEL_VISCID';
            self.data_dict.AOA=1.25;
            self.data_dict.SIDESLIP_ANGLE=0;

            self.data_dict.MACH_NUMBER=0.8;
            self.data_dict.FREESTREAM_TEMPERATURE=288.15;
            self.data_dict.FREESTREAM_PRESSURE=101325.0;

            self.data_dict.REF_ORIGIN_MOMENT_X=0.25;
            self.data_dict.REF_ORIGIN_MOMENT_Y=0.00;
            self.data_dict.REF_ORIGIN_MOMENT_Z=0.00;

            self.data_dict.ANGULAR_VELOCITY_X=0.00;
            self.data_dict.ANGULAR_VELOCITY_Y=0.00;
            self.data_dict.ANGULAR_VELOCITY_Z=0.00;

            self.data_dict.REF_LENGTH=1.00;
            self.data_dict.REF_AREA=1.00;

            self.data_dict.MARKER_ISOTHERMAL=300;
            self.data_dict.SYMMETRY='';
            self.data_dict.MARKER_MONITORING='';

            self.data_dict.MESH_FILENAME='';
            self.data_dict.MESH_FORMAT='';
            self.data_dict.MESH_SCALE=1;

            self.data_dict.OUTPUT_FILES=[];

            self.data_dict.INFORMATION=1;
        end
    end

    % parameter process funciton
    methods
        function setParameter(self,parameter,value)
            % set value of parameter
            %
            if ~isempty(value)
                if ~iscell(value), value={value}; end

                % add parameter, if is number, convert char to number
                for value_idx=1:length(value)
                    digital_value=str2double(value{value_idx});
                    if ~isnan(digital_value)
                        value{value_idx}=digital_value;
                    end
                end

                if length(value) == 1
                    self.data_dict.(parameter)=value{1};
                else
                    self.data_dict.(parameter)=value;
                end
            end
        end

        function value=getParameter(self,parameter)
            % obtain value of parameter
            %
            if isfield(self.data_dict,parameter)
                value=self.data_dict.(parameter);
            else
                value=[];
            end
        end

        function bool=isParameter(self,parameter)
            % obtain value of parameter
            %
            bool=isfield(self.data_dict,parameter);
        end
    end
end