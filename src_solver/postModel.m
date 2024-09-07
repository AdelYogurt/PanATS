function model=postModel(model)
% reconstruct element data into point data, then write into file
%
% copyright Adel 2023.03
%
config=model.config; % load config
geometry=model.geometry; % load geometry
numerics=model.numerics; % load numerics
output=model.output; % load output
INFORMATION=config.INFORMATION; % whether print information

% load geometry
pnt_list=geometry.point_list;
elem_list=geometry.element_list;
area_list=geometry.area_list;

%% average element data to point data

pnt_num=size(pnt_list,1);
elem_num=size(elem_list,1);
type_list=fieldnames(numerics);

% sum of weight
sum_weight_list=zeros(pnt_num,1);
for elem_idx=1:elem_num
    if isnan(elem_list(elem_idx,4)),node_idx_list=elem_list(elem_idx,1:3);
    else,node_idx_list=elem_list(elem_idx,1:4);end
    sum_weight_list(node_idx_list,:)=sum_weight_list(node_idx_list,:)+area_list(elem_idx);
end
sum_weight_list(sum_weight_list==0)=1;

for type_idx=1:length(type_list)
    % load element data
    type=type_list{type_idx};
    elem_data=numerics.(type);
    if size(elem_data,1) ~= elem_num,continue;end % this is not element data

    % convert element data to point data
    pnt_data=zeros(pnt_num,size(elem_data,2));
    for elem_idx=1:elem_num
        if isnan(elem_list(elem_idx,4)),node_idx_list=elem_list(elem_idx,1:3);
        else,node_idx_list=elem_list(elem_idx,1:4);end
        pnt_data(node_idx_list,:)=pnt_data(node_idx_list,:)+elem_data(elem_idx,:)*area_list(elem_idx);
    end
    pnt_data=pnt_data./sum_weight_list;

    output.(type)=pnt_data;
end

%% generate file output

%  write model into file
if isfield(config,'OUTPUT_FILES') && ~isempty(config.OUTPUT_FILES)
    out_filetype_list=config.OUTPUT_FILES;
    if ischar(out_filetype_list)
        out_filetype_list={out_filetype_list};
    end

    pnt_idx_list=(1:pnt_num)';
    % construct output data and head
    type=["""PointID""","""x""","""y""","""z"""];
    data=[pnt_idx_list,pnt_list(:,1),pnt_list(:,2),pnt_list(:,3)];

    if strcmp(model.config.SOLVER,'PANEL_INVISCID') || strcmp(model.config.SOLVER,'PANEL_VISCID')
        P_list=output.P_list;
        Cp_list=output.Cp_list;

        type=[type,"""Pressure""","""Pressure_Coefficient"""];
        data=[data,P_list,Cp_list];

        % VISCID result
        if strcmp(model.config.SOLVER,'PANEL_VISCID')
            pnt_flow_list=output.element_flow_list;
            sle_len_list=output.streamline_len_list;
            Cf_list=output.Cf_list;
            HF_list=output.HF_list;

            type=[type,"""Velocity_x""","""Velocity_y""","""Velocity_z""","""Streamline_Length"""];
            data=[data,pnt_flow_list(:,1),pnt_flow_list(:,2),pnt_flow_list(:,3),sle_len_list];
            type=[type,"""Skin_Friction_Coefficient"""];
            data=[data,Cf_list];
            type=[type,"""Heat_Flux"""];
            data=[data,HF_list];
        end
    end

    for file_idx=1:length(out_filetype_list)
        out_filetype=out_filetype_list{file_idx};

        if strcmp(out_filetype,'SURFACE_CSV')
            out_filename='surface_flow.csv';
            tbl=array2table(data,'VariableNames',type);
            writetable(tbl,out_filename)
        elseif strcmp(out_filetype,'SURFACE_TECPLOT')
            % process element list into FETRIANGLE or FEQUADRILATERAL mat
            bool=isnan(elem_list(:,4));
            if all(bool),elem_list=elem_list(:,1:3);allTri=true;
            else,elem_list(bool,4)=elem_list(bool,3);allTri=false;end

            % write tecplot data
            out_filename='surface_flow.plt';
            out_file=fopen(out_filename,'w');

            % title
            fprintf(out_file,'TITLE = "Visualization of the solution"\n');

            % variable name
            out_format=['VARIABLES =',repmat('%s,',1,length(type)-2),'%s\n'];
            fprintf(out_file,out_format,type(2:end));

            % data information
            if allTri
                fprintf(out_file,'ZONE NODES= %d, ELEMENTS= %d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n',pnt_num,elem_num);
            else
                fprintf(out_file,'ZONE NODES= %d, ELEMENTS= %d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n',pnt_num,elem_num);
            end

            % write point data
            writeMatrix(data(:,2:end),out_file);

            % write element idx
            writeMatrix(elem_list,out_file);

            fclose(out_file);
        end
    end
end

%% sort data

if config.INFORMATION
    fprintf('postModel: post model done!\n');
end

model.output=output;
end

function writeMatrix(mat,out_file)
% write matrix into file handle
%
out_format=[repmat('%.15e ',1,size(mat,2)-1),'%.15e\n'];
fprintf(out_file,out_format,mat');
end
