function postModel()
% convert face data into vertex data
%
global user_model

geometry=user_model.geometry;

dimension=geometry.dimension;
point_list=geometry.point_list;
element_list=user_model.element_list;

% load data form user_model
output_inviscid=user_model.output_inviscid;
output_streamline=user_model.output_streamline;
output_boulay=user_model.output_boulay;
output_viscid=user_model.output_viscid;
output_heat=user_model.output_heat;
output_FEM=user_model.output_FEM;
output_post=user_model.output_post;

% repeat times of vertex
point_repeat_times=zeros(size(point_list,1),1);


% load data
if ~isempty(output_inviscid)
    elem_Cp_list=output_inviscid.Cp_list;
    elem_P_list=output_inviscid.P_list;
end

if ~isempty(output_streamline)
    elem_streamline_len_list=output_streamline.streamline_len_list;
end

if ~isempty(output_boulay)
    elem_Re_x_list=output_boulay.Re_x_list;
end

if ~isempty(output_viscid)
    elem_Cf_list=output_viscid.Cf_list;
end

if ~isempty(output_heat)
    elem_HF_list=output_heat.HF_list;
end

% inviscid
if ~isempty(output_inviscid)
    Cp_list=zeros(size(point_list,1),1);
    P_list=zeros(size(point_list,1),1);
end

% streamline
if ~isempty(output_streamline)
    streamline_len_list=zeros(size(point_list,1),1);
end

% boundary layer
if ~isempty(output_boulay)
    T_2_list=zeros(size(point_list,1),1);
    P_2_list=zeros(size(point_list,1),1);
    rho_2_list=zeros(size(point_list,1),1);
    rho_e_list=zeros(size(point_list,1),1);
    V_e_list=zeros(size(point_list,1),1);
    mu_e_list=zeros(size(point_list,1),1);
    H_e_list=zeros(size(point_list,1),1);

    rho_ref_list=zeros(size(point_list,1),1);
    mu_ref_list=zeros(size(point_list,1),1);
    H_ref_list=zeros(size(point_list,1),1);
    H_r_list=zeros(size(point_list,1),1);
    Re_x_list=zeros(size(point_list,1),1);
    Re_x_ref_list=zeros(size(point_list,1),1);
end

% viscid
if ~isempty(output_viscid)
    Cf_list=zeros(size(point_list,1),1);
end

% heat
if ~isempty(output_heat)
    HF_list=zeros(size(point_list,1),1);
end

elem_num=length(element_list);
for elem_idx=1:elem_num
    elem=element_list(elem_idx);
    Node_idx=elem.Node_idx;

    if ~isempty(output_inviscid)
        Cp_list(Node_idx,:)=Cp_list(Node_idx,:)+elem_Cp_list(elem_idx);
        P_list(Node_idx,:)=P_list(Node_idx,:)+elem_P_list(elem_idx);
    end

    if ~isempty(output_streamline)
        streamline_len_list(Node_idx,:)=streamline_len_list(Node_idx,:)+elem_streamline_len_list(elem_idx);
    end

    if  ~isempty(output_boulay)
        Re_x_list(Node_idx,:)=Re_x_list(Node_idx,:)+elem_Re_x_list(elem_idx);
    end

    if ~isempty(output_viscid)
        Cf_list(Node_idx,:)=Cf_list(Node_idx,:)+elem_Cf_list(elem_idx);
    end

    if ~isempty(output_heat)
        HF_list(Node_idx,:)=HF_list(Node_idx,:)+elem_HF_list(elem_idx);
    end

    point_repeat_times(Node_idx)=point_repeat_times(Node_idx)+1;
end

elem_idx_repeat=find(point_repeat_times ~= 0);

% average velocity to each point
Cp_list(elem_idx_repeat,:)=Cp_list(elem_idx_repeat,:)./point_repeat_times(elem_idx_repeat);
P_list(elem_idx_repeat,:)=P_list(elem_idx_repeat,:)./point_repeat_times(elem_idx_repeat);
streamline_len_list(elem_idx_repeat,:)=streamline_len_list(elem_idx_repeat,:)./point_repeat_times(elem_idx_repeat);
Re_x_list(elem_idx_repeat,:)=Re_x_list(elem_idx_repeat,:)./point_repeat_times(elem_idx_repeat);
Cf_list(elem_idx_repeat,:)=Cf_list(elem_idx_repeat,:)./point_repeat_times(elem_idx_repeat);
HF_list(elem_idx_repeat,:)=HF_list(elem_idx_repeat,:)./point_repeat_times(elem_idx_repeat);

% updata data into output
output_post.Cp_list=Cp_list;
output_post.P_list=P_list;
output_post.streamline_len_list=streamline_len_list;
output_post.Re_x_list=Re_x_list;
output_post.HF_list=HF_list;
output_post.Cf_list=Cf_list;
output_post.vertex_repeat_times=point_repeat_times;

user_model.output_post=output_post;

if user_model.INFORMATION
    fprintf('postModel: post model done!\n');
end

end

