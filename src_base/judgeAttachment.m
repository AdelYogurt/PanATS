function [cross_flag,stag_flag]=judgeAttachment(nmvctr,node_list,flow_list)
% judge if this element is stagnation or have attachment line across
% 
% reference:
% [1] KENWRIGHT D N. Automatic detection of open and closed separation and
% attachment lines; proceedings of the Proceedings Visualization '98 (Cat
% No98CB36276), F 18-23 Oct. 1998, 1998 [C].
%
prec_tol=100*eps; % precision torlance
topo_tol=1e-2;

cross_flag=false(1);
stag_flag=false(1);

% project to 2D
e3=nmvctr;
e1=node_list(2,1:3)-node_list(1,1:3);
e1=e1/norm(e1);
e2=cross(e3,e1);
rot_mat=[e1;e2;e3]; % rotate matrix from global to local coordination

% transform local cooridinate to 2D cooridinate
% project flow to plane, and scale norm of flow to initial flow
flow_len_list_init=sqrt(sum(flow_list.^2,2));
flow_list=flow_list*rot_mat';
node_list=(node_list-node_list(1,1:3))*rot_mat';
flow_list(:,3)=0; % remove w velocity component
flow_len_list_proj=sqrt(sum(flow_list.^2,2));
flow_len_list_proj(flow_len_list_proj == 0)=1;
flow_list=flow_list./flow_len_list_proj.*flow_len_list_init; % extend to origin velocity length
node_list(:,end)=1;

% interpolation fit flow
coef=(node_list\flow_list(:,[1,2]))';

jacob=coef(:,1:2); % [b1,c1;b2,c2]
bias=coef(:,3); % [a1;a2]

% jacobian matrix is negative, phase is spiral
det_jacob=(jacob(1,1)*jacob(2,2)-jacob(1,2)*jacob(2,1));
if (((jacob(1,1)+jacob(2,2))^2-4*det_jacob) < -prec_tol)
    % check if spiral is out flow
    eig_val=eig(jacob);
    if real(eig_val(1)) > 0
        % out flow
        stag_pnt=-jacob\bias;
        node_proj_list=(node_list(:,1:2)-stag_pnt');
        offset=norm(max(node_proj_list)-min(node_proj_list))*topo_tol;
        node_proj_list=geomCurveOffset(node_proj_list,offset);

        % if (0,0) inside element
        if inpolygon(0,0,node_proj_list(:,1),node_proj_list(:,2))
            cross_flag=true(1);
            stag_flag=true(1);
        end
    end
    return;
end

% project element to normalized coordinates
[eig_vctr,eig_val]=eig(jacob);
eig_val=[eig_val(1,1),eig_val(2,2)];
if (eig_val(1) < eig_val(2))
    eig_val=fliplr(eig_val);
    eig_vctr=fliplr(eig_vctr);
end
% from now, eig_value(1) > eig_value(2)
% X corresponds to the eigenvectors 1
% Y corresponds to the eigenvectors 2

% evaluate eigenvalues of jacobian matrix
% if all of eigenvalues is zero, process stop
if all(abs(eig_val) < prec_tol)
    stag_pnt=[0;0];
    return;
end

% sometimes jacobian matrix can be singularity and need process respectivity
if abs(eig_val(2)) < prec_tol
    % check if cross separation line
    if abs(jacob(1)) < prec_tol
        x_c=-bias(1);
    else
        x_c=-bias(1)/jacob(1);
    end
    if abs(jacob(4)) < prec_tol
        y_c=-bias(2);
    else
        y_c=-bias(2)/jacob(4);
    end
    stag_pnt=[x_c;y_c];
elseif abs(eig_val(1)) < prec_tol
    % converge line
    if abs(jacob(1)) < prec_tol
        x_c=-bias(1);
    else
        x_c=-bias(1)/jacob(1);
    end
    if abs(jacob(4)) < prec_tol
        y_c=-bias(2);
    else
        y_c=-bias(2)/jacob(4);
    end
    stag_pnt=[x_c;y_c];
    return;
else
    if rcond(jacob) < prec_tol
        stag_pnt=lsqminnorm(jacob,-bias);
    else
        stag_pnt=-jacob\bias;
    end
end

% normalize point_2D
node_proj_list=(node_list(:,1:2)-stag_pnt')/eig_vctr';
flow_proj_list=flow_list(:,1:2)/eig_vctr';

% velocity do not too close
min_bou=min(node_proj_list,[],1);
max_bou=max(node_proj_list,[],1);
cntr_pnt=(min_bou+max_bou)/2;

% del_bou=[sum(max_bou-min_bou)/2,sum(max_bou-min_bou)/2]/2;
del_bou=(max_bou-min_bou)/2;

pnt_1=cntr_pnt+[del_bou(1),del_bou(2)];
pnt_2=cntr_pnt+[-del_bou(1),del_bou(2)];
Ve_1=eig_val.*pnt_1;
Ve_2=eig_val.*pnt_2;
cos_angle=Ve_1*Ve_2'/norm(Ve_1)/norm(Ve_2);
if cos_angle > 0.99
    return;
end

bou_node=node_proj_list(1:end,:);
if det(eig_vctr) < 0 % if det(eig_vector) less than 0, rotation
    bou_node=flipud(bou_node);
end
offset=norm(max(bou_node)-min(bou_node))*topo_tol;
bou_node=geomCurveOffset(bou_node,offset);

% judge phase portrait
% eigenvalue less than zero is concentrate
% X corresponds to the eigenvectors 1
% Y corresponds to the eigenvectors 2
if ((eig_val(1) > 0) && (eig_val(2) >= 0))
    if abs(eig_val(1)-eig_val(2)) < topo_tol
        % proper node
        % if (0,0) inside element
        if inpolygon(0,0,bou_node(:,1),bou_node(:,2))
            cross_flag=true(1);
            stag_flag=true(1);
        end
    else
        % repelling node, check small eigenvalue corresponded axis
        % if cross Y
        if judgeCrossY(bou_node(1:end,[1,2]))
            cross_flag=true(1);
        end
    end
elseif ((eig_val(1) < -0) && (eig_val(2) < -0))
    % attracting node, check
    % concentrate line is not stagnation point
    
elseif ((eig_val(1) > 0) && (eig_val(2) < -0))
    % saddle, judge axis X and axis Y which is separation line
    % means which eigenvalue is large than zero
    % repelling node, check small eigenvalue corresponded axis
    if judgeCrossY(bou_node(1:end,[1,2]))
        cross_flag=true(1);
    end
end

end

function cross_flag=judgeCrossX(node_list)
% function to judge if curve cross X
%
cross_flag=false(1);
node_num=size(node_list,1);
for node_idx=1:node_num
    node_jdx=node_idx+1;
    if node_jdx > node_num
        node_jdx=1;
    end
    if (((node_list(node_idx,2) <= 0) && ...
            (node_list(node_jdx,2) >= 0)) || ...
            ((node_list(node_idx,2) >= 0) && ...
            (node_list(node_jdx,2) <= 0)) )
        cross_flag=true(1);
        break;
    else
        cross_flag=false(1);
    end
end
end

function cross_flag=judgeCrossY(node_list)
% function to judge if curve cross Y
%
cross_flag=false(1);
node_num=size(node_list,1);
for node_idx=1:node_num
    node_jdx=node_idx+1;
    if node_jdx > node_num
        node_jdx=1;
    end
    if (((node_list(node_idx,1) <= 0) && ...
            (node_list(node_jdx,1) >= 0)) || ...
            ((node_list(node_idx,1) >= 0) && ...
            (node_list(node_jdx,1) <= 0)) )
        cross_flag=true(1);
        break;
    else
        cross_flag=false(1);
    end
end
end
