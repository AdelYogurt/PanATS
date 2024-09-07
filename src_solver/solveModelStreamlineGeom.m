function [model,SLE_out]=solveModelStreamlineGeom(model)
% identify stagnation elemenet and attachment line element
% base on geometric streamline method calculate streamline
% calculate element reference streamline length
%
% reference:
% [1] KENWRIGHT D N. Automatic detection of open and closed separation and
% attachment lines; proceedings of the Proceedings Visualization '98 (Cat
% No98CB36276), F 18-23 Oct. 1998, 1998 [C].
%
% copyright Adel 2023.03
%
config=model.config; % load config
geometry=model.geometry; % load geometry
numerics=model.numerics; % load numerics
INFORMATION=config.INFORMATION; % whether print information

% load config
AOA=config.AOA;
AOS=config.SIDESLIP_ANGLE;
SYMMETRY=config.SYMMETRY;
ref_point=[config.REF_ORIGIN_MOMENT_X,config.REF_ORIGIN_MOMENT_Y,config.REF_ORIGIN_MOMENT_Z];
omega=[config.ANGULAR_VELOCITY_X,config.ANGULAR_VELOCITY_Y,config.ANGULAR_VELOCITY_Z]/180*pi;

% load geometry
dim=geometry.dimension;
pnt_list=geometry.point_list;
elem_list=geometry.element_list;
conn_mat=geometry.connectivity_matrix;
cntr_pnt_list=geometry.center_point_list;
E_nmvctr_list=geometry.EN_vector_list;
P_nmvctr_list=geometry.PN_vector_list;
area_list=geometry.area_list;

geom_tol=1e-6;
DRAW_FIGURE_FLAG=0; % debug mode, whether draw data

% initialize numerics
numerics.element_flow_list=[];
numerics.point_flow_list=[];
numerics.element_attachment_list=[]; % element attachment list
numerics.boolean_attachment_list=[]; % whether element is attachment element
numerics.streamline_len_list=[]; % cross streamline length of element
numerics.streamline_index_list=[]; % cross streamline index of element
numerics.streamline_list=[]; % all streamline cross point
numerics.line_len_list=[]; % all streamline section length
numerics.cross_index_list=[]; % all streamline section element index

% initialize output
SLE_out.max_streamline_len=[];

%% pre process

% calculate free flow vector
coord_vec=coordVecToOri(AOA,AOS,dim);
E_FFvctr_list=coord_vec(:,1)'-cross(repmat(omega,length(elem_list),1),cntr_pnt_list-ref_point,2);
P_FFvctr_list=coord_vec(:,1)'-cross(repmat(omega,length(pnt_list),1),pnt_list-ref_point,2);

if DRAW_FIGURE_FLAG,displayModel(model);end

% calculate surface velocity of each element and point
elem_dot_norm_vec=sum(E_nmvctr_list.*E_FFvctr_list,2);
elem_flow_list=(E_FFvctr_list-E_nmvctr_list.*elem_dot_norm_vec);
pnt_dot_norm_vec=sum(P_nmvctr_list.*P_FFvctr_list,2);
pnt_flow_list=(P_FFvctr_list-P_nmvctr_list.*pnt_dot_norm_vec);

% identify attachment element
elem_num=length(elem_list);
elem_att_list=[]; % element_index
bool_att_list=false(elem_num,1);

if dim == 2

else
    for elem_idx=1:elem_num
        if isnan(elem_list(elem_idx,4)),node_idx_list=elem_list(elem_idx,1:3);
        else,node_idx_list=elem_list(elem_idx,1:4);end
        nmvctr=E_nmvctr_list(elem_idx,:);

        node_list=pnt_list(node_idx_list,1:3);
        flow_list=pnt_flow_list(node_idx_list,:);
        node_list=[node_list;cntr_pnt_list(elem_idx,1:3)];
        flow_list=[flow_list;elem_flow_list(elem_idx,:)];

        % if any(vecnorm(node_list-[0.29,0.029,0.029],2,2) < 1e-2) && any(vecnorm(node_list-[0.279288,0.029,0.0195],2,2) < 1e-2)...
        %          && any(vecnorm(node_list-[0.29,0.029,0.01],2,2) < 1e-2)
        %     disp('?');
        % end

        if dot(nmvctr,E_FFvctr_list(elem_idx,:)) < 0 % only upwind element
            [cross_flag,stag_flag]=judgeAttachment(nmvctr,node_list,flow_list);
            if cross_flag
                if DRAW_FIGURE_FLAG,hold on;patch('Faces',[1,2,3],'Vertices',pnt_list(node_idx_list,:),'FaceColor','green','FaceAlpha',0.5);end
                elem_att_list=[elem_att_list;elem_idx];
                bool_att_list(elem_idx)=true(1);
            end
        end
    end
end

% initialize data sort array
sle_len_list=zeros(elem_num,1);
sle_idx_list=zeros(elem_num,1);
sle_list={};
line_len_list={};
cross_idx_list={};
max_sle_len=0; % max length of streamline

if dim == 2

else
    % select stagnation element center point as start point to generate streamline
    for att_idx=1:length(elem_att_list)
        elem_idx=elem_att_list(att_idx);

        if sle_idx_list(elem_idx),continue;end % element with streamline donnot need search

        [streamline,line_len,cross_idx]=geomStreamline(pnt_list,elem_list,conn_mat,....
            elem_flow_list,E_nmvctr_list,elem_idx,bool_att_list);
        if DRAW_FIGURE_FLAG,line(streamline(:,1),streamline(:,2),streamline(:,3),'Color','r');end

        sle_idx=length(sle_list)+1;
        sle_list{sle_idx}=streamline;
        line_len_list{sle_idx}=line_len;
        cross_idx_list{sle_idx}=cross_idx;

        for data_idx=2:size(streamline,1)
            cross_elem_idx=cross_idx(data_idx,1);
            if sle_idx_list(cross_elem_idx)
                % if streamline length larger than exist streamline, replace it
                if sle_len_list(cross_elem_idx) > line_len(data_idx)
                    sle_len_list(cross_elem_idx)=line_len(data_idx);
                    sle_idx_list(cross_elem_idx)=sle_idx;
                end
            else
                sle_len_list(cross_elem_idx)=line_len(data_idx);
                sle_idx_list(cross_elem_idx)=sle_idx;
            end
        end

        if line_len(end) > max_sle_len
            max_sle_len=line_len(end);
        end
    end

    % adopt greedy algorithm, searching from most from free flow
    proj_cntr=sum(cntr_pnt_list.*coord_vec(:,1)',2);
    [~,order_elem_idx_list]=sort(proj_cntr);

    % check which element do not cross by streamline
    % if exist, add streamline that cross this element
    % if not, calculate new streamline that cross this element
    for search_idx=1:elem_num
        elem_idx=order_elem_idx_list(search_idx);

        if sle_idx_list(elem_idx) || abs(1-elem_dot_norm_vec(elem_idx,:)) < geom_tol
            continue; % element with streamline or have no flow donnot need search
        end

        % if donot have streamline cross, add new streamline
        [streamline,line_len,cross_idx]=geomStreamline(pnt_list,elem_list,conn_mat,....
            elem_flow_list,E_nmvctr_list,elem_idx,bool_att_list);
        if DRAW_FIGURE_FLAG,line(streamline(:,1),streamline(:,2),streamline(:,3),'Color','r');end

        sle_idx=length(sle_list)+1;
        sle_list{sle_idx}=streamline;
        line_len_list{sle_idx}=line_len;
        cross_idx_list{sle_idx}=cross_idx;

        for data_idx=2:size(streamline,1)
            cross_elem_idx=cross_idx(data_idx,1);
            if sle_idx_list(cross_elem_idx)
                % if streamline length larger than exist streamline, replace it
                if sle_len_list(cross_elem_idx) > line_len(data_idx)
                    sle_len_list(cross_elem_idx)=line_len(data_idx);
                    sle_idx_list(cross_elem_idx)=sle_idx;
                end
            else
                sle_len_list(cross_elem_idx)=line_len(data_idx);
                sle_idx_list(cross_elem_idx)=sle_idx;
            end
        end

        if line_len(end) > max_sle_len
            max_sle_len=line_len(end);
        end
    end
end

% dead water zone set to max streamline len
sle_len_list(sle_len_list == 0)=max_sle_len;

%% sort data

if config.INFORMATION
    fprintf('solveModelStreamline: streamline solve done!\n');
    fprintf('solveModelStreamline: result\n');
    fprintf('max streamline length: %14f\n',max_sle_len);
end

numerics.element_flow_list=elem_flow_list;
numerics.point_flow_list=pnt_flow_list;
numerics.element_attachment_list=elem_att_list;
numerics.boolean_attachment_list=bool_att_list;
numerics.streamline_len_list=sle_len_list;
numerics.streamline_index_list=sle_idx_list;
numerics.streamline_list=sle_list;
numerics.line_len_list=line_len_list;
numerics.cross_index_list=cross_idx_list;

SLE_out.max_streamline_len=max_sle_len;

model.numerics=numerics;

end
