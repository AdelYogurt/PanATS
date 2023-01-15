function solveModelStreamline...
    (rou_1,V_1,T_1,P_1,Ma_1,gama,AOA,SIDESLIP_ANGLE)
% calculate element reference streaamline length
%
global g_geometry g_Point g_Element g_Marker...
    ADtree_marker_element HATS_element_list...
    streamline_output

vector_flow=g_geometry.vector_flow;

geometry_torlance=1e-6;

marker_element_number=size(HATS_element_list,1);

% calculate streamline length
inflow_direction_list=zeros(marker_element_number,3);
streamline_length_list=zeros(marker_element_number,1);

bottom_index_list=[];
for element_index=1:marker_element_number
    center_point=ADtree_marker_element.center_point_list(element_index,:);
    out_index_list=HATS_element_list(element_index).out_index_list;
    normal_vector=HATS_element_list(element_index).normal_vector;
    
    if abs(normal_vector*vector_flow-1) < geometry_torlance
        bottom_index_list=[bottom_index_list,element_index];
        continue;
    end
    
    inflow_direction=inflow_direction_list(element_index,:);
    streamline_length=streamline_length_list(element_index);
    
    for out_node_index__=1:length(out_index_list)
        out_node_index=out_index_list(out_node_index__);
        
        out_center_point=ADtree_marker_element.center_point_list(out_node_index,:);
        
        % base on normal_vector calculate inflow_direction and streamline_length
        inflow_direction_out=(out_center_point-center_point);
        streamline_length_out=sqrt(sum(inflow_direction_out.^2));
        inflow_direction_out=inflow_direction_out/streamline_length_out;
        
        % correct length and vector base on parent inflow_direction and streamline_length
        if ~isempty(inflow_direction)
            inflow_direction_out=...
                inflow_direction_out*streamline_length_out+...
                inflow_direction*streamline_length; % vector plus
            streamline_length_out=sqrt(sum(inflow_direction_out.^2));
            inflow_direction_out=inflow_direction_out/streamline_length_out;
        end
        
        % compare exist steamline, if shorter than repace
        if (streamline_length_list(out_node_index) ~= 0)
            if streamline_length_list(out_node_index) > streamline_length_out
                inflow_direction_list(out_node_index,:)=inflow_direction_out;
                streamline_length_list(out_node_index)=streamline_length_out;
            end
        else
            inflow_direction_list(out_node_index,:)=inflow_direction_out;
            streamline_length_list(out_node_index)=streamline_length_out;
        end
    end
end
inflow_direction_list(bottom_index_list,:)=repmat(vector_flow',length(bottom_index_list),1);
streamline_length_list(bottom_index_list)=max(streamline_length_list);

streamline_output.inflow_direction_list=inflow_direction_list;
streamline_output.streamline_length_list=streamline_length_list;
end