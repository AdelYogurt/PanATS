function postModel()
% reconstruct cell data into vertex data, then write into file
%
% copyright Adel 2023.03
%
global user_model

config=user_model.config;
geometry=user_model.geometry;
element_list=user_model.element_list;
SYMMETRY=config.SYMMETRY;

% load geometry
dimension=geometry.dimension;
point_list=geometry.point_list;
area_list=geometry.area_list;

% load data form user_model
output_inviscid=user_model.output_inviscid;
output_streamline=user_model.output_streamline;
output_boulay=user_model.output_boulay;
output_viscid=user_model.output_viscid;
output_heat=user_model.output_heat;

output_post=user_model.output_post;

pnt_num=size(point_list,1);
% inviscid
if ~isempty(output_inviscid)
    elem_Cp_list=output_inviscid.Cp_list;
    elem_P_list=output_inviscid.P_list;
    Cp_list=zeros(pnt_num,1);
    P_list=zeros(pnt_num,1);
end

% streamline
if ~isempty(output_streamline)
    elem_streamline_len_list=output_streamline.streamline_len_list;
    elem_surface_flow_list=output_streamline.surface_flow_list;
    streamline_len_list=zeros(pnt_num,1);
    surface_flow_list=zeros(pnt_num,3);
end

% viscid
if ~isempty(output_viscid)
    elem_Cf_list=output_viscid.Cf_list;
    Cf_list=zeros(pnt_num,1);
end

% heat
if ~isempty(output_heat)
    elem_HF_list=output_heat.HF_list;
    HF_list=zeros(pnt_num,1);
end

area_sum_list=zeros(pnt_num,1); % repeat times of area sum
elem_num=length(element_list);
for elem_idx=1:elem_num
    elem=element_list(elem_idx);
    Node_idx=elem.Node_idx;
    area=area_list(elem_idx);

    % add data to each point
    if ~isempty(output_inviscid)
        Cp_list(Node_idx,:)=Cp_list(Node_idx,:)+elem_Cp_list(elem_idx,:)*area;
        P_list(Node_idx,:)=P_list(Node_idx,:)+elem_P_list(elem_idx,:)*area;
    end

    if ~isempty(output_streamline)
        streamline_len_list(Node_idx,:)=streamline_len_list(Node_idx,:)+elem_streamline_len_list(elem_idx,:)*area;
        surface_flow_list(Node_idx,:)=surface_flow_list(Node_idx,:)+elem_surface_flow_list(elem_idx,:)*area;
    end

    if ~isempty(output_viscid)
        Cf_list(Node_idx,:)=Cf_list(Node_idx,:)+elem_Cf_list(elem_idx,:)*area;
    end

    if ~isempty(output_heat)
        HF_list(Node_idx,:)=HF_list(Node_idx,:)+elem_HF_list(elem_idx,:)*area;
    end

    area_sum_list(Node_idx)=area_sum_list(Node_idx)+area;
end

area_sum_list(area_sum_list == 0)=1;

% average velocity to each point and updata data into output
if ~isempty(output_inviscid)
    Cp_list=Cp_list./area_sum_list;
    P_list=P_list./area_sum_list;
    output_post.Cp_list=Cp_list;
    output_post.P_list=P_list;
end

if ~isempty(output_streamline)
    streamline_len_list=streamline_len_list./area_sum_list;
    surface_flow_list=surface_flow_list./area_sum_list;
    output_post.streamline_len_list=streamline_len_list;
    output_post.surface_flow_list=surface_flow_list;
end

if ~isempty(output_viscid)
    Cf_list=Cf_list./area_sum_list;
    output_post.Cf_list=Cf_list;
end

if ~isempty(output_heat)
    HF_list=HF_list./area_sum_list;
    output_post.HF_list=HF_list;
end

user_model.output_post=output_post;

% write model into file
if isfield(config,'OUTPUT_FILES') && ~isempty(config.OUTPUT_FILES)
    out_filetype_list=config.OUTPUT_FILES;
    if ischar(out_filetype_list)
        out_filetype_list={out_filetype_list};
    end

    pnt_idx_list=(1:pnt_num)';
    % construct output data and head
    type=["""PointID""","""x""","""y""","""z"""];
    data=[pnt_idx_list,point_list(:,1),point_list(:,2),point_list(:,3)];

    % inviscid
    if ~isempty(output_inviscid)
        type=[type,"""Pressure""","""Pressure_Coefficient"""];
        data=[data,P_list,Cp_list];
    end

    % streamline
    if ~isempty(output_streamline)
        type=[type,"""Velocity_x""","""Velocity_y""","""Velocity_z""","""Streamline_Length"""];
        data=[data,surface_flow_list(:,1),surface_flow_list(:,2),surface_flow_list(:,3),streamline_len_list];
    end

    % viscid
    if ~isempty(output_viscid)
        type=[type,"""Skin_Friction_Coefficient"""];
        data=[data,Cf_list];
    end

    % heat
    if ~isempty(output_heat)
        type=[type,"""Heat_Flux"""];
        data=[data,HF_list];
    end

    for file_idx=1:length(out_filetype_list)
        out_filetype=out_filetype_list{file_idx};

        if strcmp(out_filetype,'SURFACE_CSV')
            out_filename='surface_flow.csv';
            tbl=array2table(data,'VariableNames',type);
            writetable(tbl,out_filename)
        elseif strcmp(out_filetype,'SURFACE_TECPLOT')
            % process element index into FETRIANGLE or FEQUADRILATERAL mat
            elem_num=length(element_list);
            each_node_num=3;
            elem_idx_mat=zeros(elem_num,each_node_num);
            allTri=true;
            for elem_idx=1:elem_num
                elem=element_list(elem_idx);
                Node_idx=elem.Node_idx;

                if allTri && each_node_num < elem.node_num
                    % change tri mat into qual mat
                    elem_idx_mat=[elem_idx_mat,zeros(elem_num,1)];
                    elem_idx_mat(1:elem_idx,end)=elem_idx_mat(1:elem_idx,end-1);
                    each_node_num=4;
                    allTri=false;
                end

                if elem.node_num < each_node_num, Node_idx=[Node_idx,Node_idx(end)];end

                elem_idx_mat(elem_idx,:)=Node_idx(1:each_node_num);
            end

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
                fprintf(out_file,'ZONE NODES= %d, ELEMENTS= %d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE',pnt_num,elem_num);
            else
                fprintf(out_file,'ZONE NODES= %d, ELEMENTS= %d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL',pnt_num,elem_num);
            end

            fclose(out_file);

            % write point data
            writematrix(data(:,2:end),out_filename,'WriteMode','append','FileType','text','Delimiter',' ');

            % write element idx
            writematrix(elem_idx_mat,out_filename,'WriteMode','append','FileType','text','Delimiter',' ');
        end
    end
end

if config.INFORMATION
    fprintf('postModel: post model done!\n');
end

end

