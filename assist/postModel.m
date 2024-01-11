function postModel()
% convert face data into vertex data
%
global user_model

dimension=user_model.geometry.dimension;

point_list=user_model.point_list;
marker_list=user_model.marker_list;

MARKER_MONITERING=user_model.MARKER_MONITORING;

% load data form user_model
output_inviscid=user_model.output_inviscid;
output_streamline=user_model.output_streamline;
output_boulay=user_model.output_boulay;
output_viscid=user_model.output_viscid;
output_heat=user_model.output_heat;

output_post=user_model.output_post;

% repeat times of vertex
if ~isfield(output_post,'vertex_repeat_times')
    vertex_flag=true(1);
    vertex_repeat_times=[zeros(size(point_list,1),1)];
else
    vertex_flag=false(1);
    vertex_repeat_times=output_post.vertex_repeat_times;
    marker_boolean=output_post.marker_boolean;
end

% load data
if ~isempty(output_inviscid)
    marker_Cp_list=output_inviscid.Cp_list;
    marker_P_list=output_inviscid.P_list;
    marker_dFn_list=output_inviscid.dFn_list;
    marker_dMn_list=output_inviscid.dMn_list;
end

if ~isempty(output_streamline)
    marker_streamline_len_list=output_streamline.streamline_len_list;
end

if ~isempty(output_boulay)
    marker_Re_x_list=output_boulay.Re_x_list;
end

if ~isempty(output_viscid)
    marker_Cf_list=output_viscid.Cf_list;
    marker_dFs_list=output_viscid.dFs_list;
    marker_dMs_list=output_viscid.dMs_list;
end

if ~isempty(output_heat)
    marker_HF_list=output_heat.HF_list;
end

% initialize data array
total_point_number=size(point_list,1);

Cp_list=zeros(total_point_number,1);
P_list=zeros(total_point_number,1);
dFn_list=zeros(total_point_number,1);
dMn_list=zeros(total_point_number,1);

streamline_len_list=zeros(total_point_number,1);

HF_list=zeros(total_point_number,1);
Re_x_list=zeros(total_point_number,1);

Cf_list=zeros(total_point_number,1);
dFs_list=zeros(total_point_number,1);
dMs_list=zeros(total_point_number,1);

for moniter_index=1:length(MARKER_MONITERING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITERING{moniter_index},marker_list);
    for element_index=1:marker_list(marker_index).element_number
        element=marker_element(element_index);

        point_index_list=element.point_index_list;

        if ~isempty(output_inviscid)
            Cp_list(point_index_list,:)=Cp_list(point_index_list,:)+marker_Cp_list{marker_index}(element_index);
            P_list(point_index_list,:)=P_list(point_index_list,:)+marker_P_list{marker_index}(element_index);
            dFn_list(point_index_list,:)=dFn_list(point_index_list,:)+marker_dFn_list{marker_index}(element_index);
            dMn_list(point_index_list,:)=dMn_list(point_index_list,:)+marker_dMn_list{marker_index}(element_index);
        end

        if ~isempty(output_streamline)
            streamline_len_list(point_index_list,:)=streamline_len_list(point_index_list,:)+marker_streamline_len_list{marker_index}(element_index);
        end

        if ~isempty(output_heat)
            HF_list(point_index_list,:)=HF_list(point_index_list,:)+marker_HF_list{marker_index}(element_index);
            Re_x_list(point_index_list,:)=Re_x_list(point_index_list,:)+marker_Re_x_list{marker_index}(element_index);
        end

        if ~isempty(output_viscid)
            Cf_list(point_index_list,:)=Cf_list(point_index_list,:)+marker_Cf_list{marker_index}(element_index);
            dFs_list(point_index_list,:)=dFs_list(point_index_list,:)+marker_dFs_list{marker_index}(element_index);
            dMs_list(point_index_list,:)=dMs_list(point_index_list,:)+marker_dMs_list{marker_index}(element_index);
        end

        if vertex_flag
            vertex_repeat_times(point_index_list)=vertex_repeat_times(point_index_list)+1;
        end
    end
end
if vertex_flag
    marker_boolean=vertex_repeat_times ~= 0;
end
% average velocity to each point
Cp_list(marker_boolean,:)=Cp_list(marker_boolean,:)./vertex_repeat_times(marker_boolean);
P_list(marker_boolean,:)=P_list(marker_boolean,:)./vertex_repeat_times(marker_boolean);
dFn_list(marker_boolean,:)=dFn_list(marker_boolean,:)./vertex_repeat_times(marker_boolean);
dMn_list(marker_boolean,:)=dMn_list(marker_boolean,:)./vertex_repeat_times(marker_boolean);

streamline_len_list(marker_boolean,:)=streamline_len_list(marker_boolean,:)./vertex_repeat_times(marker_boolean);

HF_list(marker_boolean,:)=HF_list(marker_boolean,:)./vertex_repeat_times(marker_boolean);
Re_x_list(marker_boolean,:)=Re_x_list(marker_boolean,:)./vertex_repeat_times(marker_boolean);

Cf_list(marker_boolean,:)=Cf_list(marker_boolean,:)./vertex_repeat_times(marker_boolean);
dFs_list(marker_boolean,:)=dFs_list(marker_boolean,:)./vertex_repeat_times(marker_boolean);
dMs_list(marker_boolean,:)=dMs_list(marker_boolean,:)./vertex_repeat_times(marker_boolean);

% updata data into output
output_post.Cp_list=Cp_list;
output_post.P_list=P_list;
output_post.dFn_list=dFn_list;
output_post.dMn_list=dMn_list;

output_post.streamline_len_list=streamline_len_list;

output_post.HF_list=HF_list;
output_post.Re_x_list=Re_x_list;

output_post.Cf_list=Cf_list;
output_post.dFs_list=dFs_list;
output_post.dMs_list=dMs_list;

output_post.dMs_list=dMs_list;

if vertex_flag
    output_post.vertex_repeat_times=vertex_repeat_times;
    output_post.marker_boolean=marker_boolean;
end

user_model.output_post=output_post;

if user_model.INFORMATION
    fprintf('postModel: post model done!\n');
end

end

