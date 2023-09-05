function [max_streamline_len]=solveModelStreamline()
% identify specifical element(stagnation elemenet, separate element)
% specifical element is the start of streamline
% calculate element reference streamline length
%
% reference: [1] KENWRIGHT D N. Automatic detection of open and closed
% separation and attachment lines; proceedings of the Proceedings
% Visualization '98 (Cat No98CB36276), F 18-23 Oct. 1998, 1998 [C].
%
% copyright Adel 2023.03
%
global user_model
DRAW_FIGURE_FLAG=0; % whether draw data

geometry_torlance=1e-9;

geometry=user_model.geometry;
element_list=user_model.element_list;
SYMMETRY=user_model.SYMMETRY;

dimension=geometry.dimension;
point_list=geometry.point_list;

% load geometry
center_point_list=geometry.center_point_list;
area_list=geometry.area_list;
normal_vector_list=geometry.normal_vector_list;

% calculate free flow vector
free_flow_vector=calFreeFlowDirection(user_model.AOA,user_model.SIDESLIP_ANGLE);
user_model.free_flow_vector=free_flow_vector;

switch SYMMETRY
    case 'XOY'
        flip_operate=[1,1,-1];
        identify_dimension=3;
    case 'YOZ'
        flip_operate=[-1,1,1];
        identify_dimension=1;
    case 'ZOX'
        flip_operate=[1,-1,1];
        identify_dimension=2;
    otherwise
        error('solveModelStreamline: nuknown SYMMETRY type');
end

if DRAW_FIGURE_FLAG,displayMarker();end

% elem_num=length(element_list);
% point_num=size(point_list,1);
% % calculate normal vector of point
% point_normal_vector_list=zeros(point_num,3); % surface_flow
% vertex_repeat_times=zeros(point_num,1); % repeat times of vertex
% for elem_idx=1:elem_num
%     elem=element_list(elem_idx);
%     Node_idx=elem.Node_idx;
%     area=area_list(elem_idx);
%     normal_vector=normal_vector_list(elem_idx,:);
% 
%     % add data to each point
%     point_normal_vector_list(Node_idx,:)=...
%         point_normal_vector_list(Node_idx,:)+normal_vector*area;
%     vertex_repeat_times(Node_idx)=...
%         vertex_repeat_times(Node_idx)+1;
% end
% % identify which point of element on symmetry face
% % add symmetry velocity to point
% % (velocity identify_dimension set to zero)
% Bool_sym=abs(point_list(:,identify_dimension)) <= geometry_torlance;
% point_normal_vector_list(Bool_sym,identify_dimension)=0;
% point_normal_vector_list=point_normal_vector_list./vecnorm(point_normal_vector_list,2,2);

% calculate surface velocity of each element
dot_nor_vec=normal_vector_list*free_flow_vector;
elem_flow_list=(free_flow_vector'-normal_vector_list.*dot_nor_vec);

% % calculate surface velocity of each point
% point_dot_nor_vec=point_normal_vector_list*free_flow_vector;
% point_flow_list=(free_flow_vector'-point_normal_vector_list.*point_dot_nor_vec);

elem_num=length(element_list);
point_num=size(point_list,1);
% calculate normal vector of point
point_flow_list=zeros(point_num,3); % surface_flow
area_sum_list=zeros(point_num,1); % repeat times of vertex
for elem_idx=1:elem_num
    elem=element_list(elem_idx);
    Node_idx=elem.Node_idx;
    area=area_list(elem_idx);
    elem_flow=elem_flow_list(elem_idx,:);

    if norm(elem_flow) > geometry_torlance
        % add data to each point
        point_flow_list(Node_idx,:)=...
            point_flow_list(Node_idx,:)+elem_flow*area;
        area_sum_list(Node_idx)=...
            area_sum_list(Node_idx)+area;
    end
end
% identify which point of element on symmetry face
% add symmetry velocity to point
% (velocity identify_dimension set to zero)
Bool_sym=abs(point_list(:,identify_dimension)) <= geometry_torlance;
point_flow_list(Bool_sym,identify_dimension)=0;
area_sum_list(area_sum_list == 0)=1;
point_flow_list=point_flow_list./area_sum_list;

% record
user_model.output_post.surface_flow_list=point_flow_list;

elem_num=length(element_list);
% identify attachment element
elem_attach_list=[]; % element_index
bool_attach_list=false(elem_num,1);
for elem_idx=1:elem_num
    elem=element_list(elem_idx);
    Node_idx=elem.Node_idx;
    normal_vector=normal_vector_list(elem_idx,:);

%     node_list=point_list(Node_idx,:);
%     node_flow_list=point_flow_list(Node_idx,:);
%     node_list=[node_list;center_point_list(elem_idx,:)];
%     node_flow_list=[node_flow_list;elem_flow_list(elem_idx,:)];

    Vertex_next=elem.Vertex_next;
    elem_idx_list=[elem_idx,Vertex_next];
    elem_idx_list(elem_idx_list==0)=[];
    node_list=center_point_list(elem_idx_list,:);
    node_flow_list=elem_flow_list(elem_idx_list,:);
    Bool=vecnorm(node_flow_list,2,2) < geometry_torlance;
    node_list(Bool,:)=[];
    node_flow_list(Bool,:)=[];

    if dot(normal_vector,free_flow_vector) < 0 % only upwind element
        [cross_flag,stagnation_flag]=judgeAttachment(normal_vector,node_list,node_flow_list);
        if cross_flag
            elem_attach_list=[elem_attach_list,elem_idx];
            bool_attach_list(elem_idx)=true(1);

            if DRAW_FIGURE_FLAG,hold on;patch('Faces',[1,2,3],'Vertices',point_list(Node_idx,:),'FaceColor','green','FaceAlpha',0.5);end
        end
    end
end

elem_num=length(element_list);
% initialize data sort array
% streamline_len_list, streamline_list
streamline_len_list=zeros(elem_num,1);
streamline_idx_list=zeros(elem_num,1);
streamline_list={};
line_len_list={};
cross_idx_list={};
max_streamline_len=0; % max length of streamline

% select stagnation element center point as start point to generate streamline
for attachment_idx=1:length(elem_attach_list)
    elem_idx=elem_attach_list(attachment_idx);

    if ~streamline_idx_list(elem_idx)
        streamline_idx=length(streamline_list)+1;
        [streamline,line_len,cross_idx]=calStreamline(elem_idx,...
            point_list,element_list,bool_attach_list,elem_flow_list,normal_vector_list);
        if DRAW_FIGURE_FLAG,line(streamline(:,1),streamline(:,2),streamline(:,3),'Color','r');end

        streamline_list{streamline_idx}=streamline;
        line_len_list{streamline_idx}=line_len;
        cross_idx_list{streamline_idx}=cross_idx;

        for data_idx=2:size(streamline,1)
            cross_elem_idx=cross_idx(data_idx,1);
            if streamline_idx_list(cross_elem_idx)
                % if streamline length larger than exist streamline, replace it
                if streamline_len_list(cross_elem_idx) > line_len(data_idx)
                    streamline_len_list(cross_elem_idx)=line_len(data_idx);
                    streamline_idx_list(cross_elem_idx)=streamline_idx;
                end
            else
                streamline_len_list(cross_elem_idx)=line_len(data_idx);
                streamline_idx_list(cross_elem_idx)=streamline_idx;
            end
        end

        if line_len(end) > max_streamline_len
            max_streamline_len=line_len(end);
        end
    end
end

% check which element do not cross by streamline
% if exist, add streamline that cross this element
% if not, calculate new streamline that cross this element
for elem_idx=1:elem_num
    if ~streamline_idx_list(elem_idx) && ...
            abs(1-dot_nor_vec(elem_idx,:)) >= geometry_torlance
        % if donot have streamline cross, add new streamline
        streamline_idx=length(streamline_list)+1;
        [streamline,line_len,cross_idx]=calStreamline(elem_idx,...
            point_list,element_list,bool_attach_list,elem_flow_list,normal_vector_list);
        if DRAW_FIGURE_FLAG,line(streamline(:,1),streamline(:,2),streamline(:,3),'Color','r');end

        streamline_list{streamline_idx}=streamline;
        line_len_list{streamline_idx}=line_len;
        cross_idx_list{streamline_idx}=cross_idx;

        for data_idx=2:size(streamline,1)
            cross_elem_idx=cross_idx(data_idx,1);
            if streamline_idx_list(cross_elem_idx)
                % if streamline length larger than exist streamline, replace it
                if streamline_len_list(cross_elem_idx) > line_len(data_idx)
                    streamline_len_list(cross_elem_idx)=line_len(data_idx);
                    streamline_idx_list(cross_elem_idx)=streamline_idx;
                end
            else
                streamline_len_list(cross_elem_idx)=line_len(data_idx);
                streamline_idx_list(cross_elem_idx)=streamline_idx;
            end
        end

        if line_len(end) > max_streamline_len
            max_streamline_len=line_len(end);
        end
    end
end

output_streamline.surface_flow_list=elem_flow_list;
output_streamline.element_attach_list=elem_attach_list;
output_streamline.boolean_attach_list=bool_attach_list;
output_streamline.streamline_len_list=streamline_len_list;
output_streamline.streamline_idx_list=streamline_idx_list;
output_streamline.streamline_list=streamline_list;
output_streamline.line_len_list=line_len_list;
output_streamline.cross_idx_list=cross_idx_list;

user_model.output_streamline=output_streamline;

if user_model.INFORMATION
    fprintf('solveModelStreamline: streamline solve done!\n');
    fprintf('solveModelStreamline: result\n');
    fprintf('max streamline length: %14f\n',max_streamline_len);
end

end
