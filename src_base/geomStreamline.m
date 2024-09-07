function [streamline,line_len,cross_idx]=geomStreamline(pnt_list,elem_list,conn_mat,....
    elem_flow_list,E_nmvctr_list,elem_idx,bool_att_list)
% calculate streamline from point_start and element
%
DRAW_FIGURE_FLAG=0; % debug mode, whether draw data

geom_tol=1e-6;

elem_num=length(elem_list);
max_size=elem_num;

streamline=zeros(max_size,3); % cross point
line_len=zeros(max_size,1); % current length
cross_idx=zeros(max_size,1,'uint32'); % [cross elem_idx, cross edge_idx]

if isnan(elem_list(elem_idx,4)),node_idx_list=elem_list(elem_idx,1:3);node_num=3;
else,node_idx_list=elem_list(elem_idx,1:4);node_num=4;end
node_list=pnt_list(node_idx_list,1:3);
pnt_start=sum(node_list,1)/node_num;
data_idx=1;

streamline(data_idx,:)=pnt_start;
cross_idx(data_idx,:)=elem_idx;

% if DRAW_FIGURE_FLAG,displayModel();end
if DRAW_FIGURE_FLAG,line(pnt_start(:,1),pnt_start(:,2),pnt_start(:,3),'Marker','o');end

%% upwind search begin

while elem_idx ~= 0 && data_idx <= max_size/2
    % load element information
    if isnan(elem_list(elem_idx,4)),node_idx_list=elem_list(elem_idx,1:3);node_num=3;
    else,node_idx_list=elem_list(elem_idx,1:4);node_num=4;end
    elem_flow=elem_flow_list(elem_idx,:);
    nmvctr=E_nmvctr_list(elem_idx,:);
    elem_flow_len=norm(elem_flow);
    if elem_flow_len <= geom_tol,break;end % search end in no flow element
    elem_flow=elem_flow/elem_flow_len;
    node_list=pnt_list(node_idx_list,:);

    elem_flow=-elem_flow; % notice upwind search

    % serach flow upwind end point
    [pnt_end,node_idx,node_jdx,step_len]=calElementCrossPoint...
        (pnt_start,node_list,node_num,nmvctr,elem_flow,geom_tol);
    if step_len < geom_tol % if step_len too small, try upwind by node point
        vctr_len_list=(node_list-pnt_start)*elem_flow';
        [vctr_len_max,node_max_idx]=max(vctr_len_list);node_max_idx=node_max_idx(1);
        if vctr_len_max >= geom_tol
            pnt_end=node_list(node_max_idx,:);
            step_len=norm(pnt_start-pnt_end);
        end
    end

    if step_len < geom_tol
        % cannot find upwind end point, stop
        break;
    end

    % record parameter
    data_idx=data_idx+1;
    streamline(data_idx,:)=pnt_end;
    line_len(data_idx)=line_len(data_idx-1)+step_len;
    cross_idx(data_idx)=elem_idx;
    % line(point_end(:,1),point_end(:,2),point_end(:,3),'Marker','o','Color','r')

    % end in attach element
    if bool_att_list(elem_idx)
        break;
    end
   
    % search next element
    % check if is point_start overlap node list
    node_overlap_idx=find(sum(abs(node_list-pnt_end),2) < geom_tol,1);
    if ~isempty(node_overlap_idx)
        % if overlap, we need to find next element that can upwind
        pnt_overlap_idx=node_idx_list(node_overlap_idx);
        pnt_arou_idx_list=find(conn_mat(pnt_overlap_idx,:) ~= 0);

        elem_idx_arou_list=full(conn_mat(pnt_overlap_idx,pnt_arou_idx_list));
        pnt_arou_list=pnt_list(pnt_arou_idx_list,1:3);

        vctr_arou_list=(pnt_arou_list-pnt_end);
        vctr_arou_list=vctr_arou_list./sqrt(sum(vctr_arou_list.^2,2));
        cosang_arou_list=vctr_arou_list*elem_flow';
        [~,upwind_idx]=max(cosang_arou_list);

        elem_idx=elem_idx_arou_list(upwind_idx);
    else
        % if not overlap, use node_idx to find next element to upwind
        elem_new_idx=full(conn_mat(node_idx_list(node_jdx),node_idx_list(node_idx)));
        if elem_new_idx ~= 0,elem_idx=elem_new_idx;end
    end
    pnt_start=pnt_end;

    if DRAW_FIGURE_FLAG,line(pnt_start(:,1),pnt_start(:,2),pnt_start(:,3),'Marker','o');end

    if data_idx == max_size/2
        warning('geomStreamline: upwind search cross too many element');
    end
end


%% reverse and connect point

if data_idx > 1
    % reverse streamline
    streamline(1:data_idx,:)=flipud(streamline(1:data_idx,:));
    line_len(1:data_idx)=flipud(line_len(1:data_idx));
    line_len(1:data_idx)=line_len(1)-line_len(1:data_idx);

    % reverse cross_list
    cross_idx(1:data_idx)=flipud(cross_idx(1:data_idx));
    cross_idx(2:data_idx)=cross_idx(1:data_idx-1);

    data_idx=data_idx-1;
    pnt_start=streamline(data_idx,1:3);
    elem_idx=cross_idx(data_idx+1);
end

%% downwind search begin

while elem_idx ~= 0 && data_idx <= max_size
    % load element information
    if isnan(elem_list(elem_idx,4)),node_idx_list=elem_list(elem_idx,1:3);node_num=3;
    else,node_idx_list=elem_list(elem_idx,1:4);node_num=4;end
    elem_flow=elem_flow_list(elem_idx,:);
    nmvctr=E_nmvctr_list(elem_idx,:);
    elem_flow_len=norm(elem_flow);
    if elem_flow_len <= geom_tol,break;end % search end in no flow element
    elem_flow=elem_flow/elem_flow_len;
    node_list=pnt_list(node_idx_list,:);

    % calculate next cross point
    [pnt_end,node_idx,node_jdx,step_len]=calElementCrossPoint...
        (pnt_start,node_list,node_num,nmvctr,elem_flow,geom_tol);

    if step_len < geom_tol
        % check if is point_start overlap node point cause this problem
        bool_overlap_list=sum(abs(node_list-pnt_start),2) < geom_tol;
        if any(bool_overlap_list)
            % if yes, search element arounded that can continue flow
            elem_idx_orig=elem_idx;
            node_idx=find(bool_overlap_list,1);
            pnt_idx=node_idx_list(node_idx); % overlap point index

            % loop all element around of overlap point to find if downwind
            found_flag=false;
            elem_idx_arou_list=full(conn_mat(pnt_idx,conn_mat(pnt_idx,:) ~= 0));
            for loop_idx=1:length(elem_idx_arou_list)
                elem_idx=elem_idx_arou_list(loop_idx);
                if elem_idx_orig == elem_idx,continue;end % jump self element

                % load element information
                if isnan(elem_list(elem_idx,4)),node_idx_list=elem_list(elem_idx,1:3);node_num=3;
                else,node_idx_list=elem_list(elem_idx,1:4);node_num=4;end
                elem_flow=elem_flow_list(elem_idx,:);
                nmvctr=E_nmvctr_list(elem_idx,:);
                elem_flow_len=norm(elem_flow);
                if elem_flow_len <= geom_tol,break;end % search end in no flow element
                elem_flow=elem_flow/elem_flow_len;
                node_list=pnt_list(node_idx_list,:);

                [pnt_end,node_idx,node_jdx,step_len]=calElementCrossPoint...
                    (pnt_start,node_list,node_num,nmvctr,elem_flow,geom_tol);

                if step_len >= geom_tol
                    % next element found, contiune searching
                    found_flag=true;
                    break;
                end
            end

            if ~found_flag
                % if all around element searched but there is not element found
                break;
            end
        else
            % if no, means encounter the opposite flow, stop search
            break;
        end
    end

    % record parameter
    data_idx=data_idx+1;
    streamline(data_idx,:)=pnt_end;
    line_len(data_idx)=line_len(data_idx-1)+step_len;
    cross_idx(data_idx)=elem_idx;
    %     line(point_end(:,1),point_end(:,2),point_end(:,3),'Marker','o','Color','g')

    % search next element
    elem_idx=full(conn_mat(node_idx_list(node_jdx),node_idx_list(node_idx)));
    pnt_start=pnt_end;

    if DRAW_FIGURE_FLAG,line(pnt_start(:,1),pnt_start(:,2),pnt_start(:,3),'Marker','o');end

    if data_idx == max_size
        warning('geomStreamline: downwind search cross too many element');
    end
end

streamline=streamline(1:data_idx,:);
line_len=line_len(1:data_idx,:);
cross_idx=cross_idx(1:data_idx);
end

function [pnt_end,node_idx,node_jdx,d_len]=calElementCrossPoint...
    (pnt_start,node_list,node_num,nmvctr,elem_flow,geom_tol)
% calculate next cross point
%

% elem_flow as x, normal_vector as z, point_start is origin, project point
y_flow=cross(nmvctr,elem_flow);
rot_mat=[elem_flow;y_flow;nmvctr]';
proj_node_list=(node_list-pnt_start)*rot_mat;

% base on cross x axis to calculate projected point_end
pnt_end=[min(proj_node_list(:,1)),0,0];
node_idx=0;node_jdx=0;
for node_slf_idx=1:node_num
    if node_slf_idx == node_num
        node_ref_idx=1;
    else
        node_ref_idx=node_slf_idx+1;
    end

    % calculate each line intersection of x axis
    vtx_slf_x=proj_node_list(node_slf_idx,1);
    vtx_slf_y=proj_node_list(node_slf_idx,2);
    vtx_ref_x=proj_node_list(node_ref_idx,1);
    vtx_ref_y=proj_node_list(node_ref_idx,2);

    if (vtx_slf_y <= geom_tol && vtx_ref_y >= -geom_tol) ||...
            (vtx_slf_y >= -geom_tol && vtx_ref_y <= geom_tol)
        % edge should cross positive direction of axis x
        if abs(vtx_slf_y-vtx_ref_y) <= geom_tol
            % if this edge lay on axis x, then judge vtx_ref_x
            if vtx_ref_x > pnt_end(1)
                pnt_end=[vtx_ref_x,0,proj_node_list(node_ref_idx,3)];
                node_idx=node_slf_idx;
                node_jdx=node_ref_idx;
            end
        else
            cross_x=vtx_slf_x-vtx_slf_y*(vtx_slf_x-vtx_ref_x)/(vtx_slf_y-vtx_ref_y);

            if size(node_list,1) == 4
                % linear interp coordination of corss y
                pct=(cross_x-vtx_slf_x)/(vtx_ref_x-vtx_slf_x);
                cross_z=(1-pct)*proj_node_list(node_slf_idx,3)+pct*proj_node_list(node_ref_idx,3);
            else
                cross_z=0;
            end

            if cross_x > pnt_end(1)
                pnt_end=[cross_x,0,cross_z];
                node_idx=node_slf_idx;
                node_jdx=node_ref_idx;
            end
        end
    end
end

if node_idx
    % re project to real coordinate
    pnt_end=pnt_end*rot_mat'+pnt_start;
    d_len=norm(pnt_start-pnt_end);
else
    d_len=0;
    pnt_end=[];
end
end

