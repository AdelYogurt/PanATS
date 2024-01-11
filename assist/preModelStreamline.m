function preModelStreamline()
% base on half edge data structure to distribute nearby element of element
% method is basing on vertex.neaby_element_list disturbuting nearby element
%
% copyright Adel 2023.03
%
global user_model

edge_empty=user_model.edge_empty;
edge_list=user_model.edge_list;

marker_list=user_model.marker_list;

MARKER_MONITERING=user_model.MARKER_MONITORING;

% calculate free flow vector
free_flow_vector=[1;0;0];

AOA=user_model.AOA/180*pi;
cos_AOA=cos(AOA);
sin_AOA=sin(AOA);
rotation_AOA=[
    cos_AOA 0 -sin_AOA;
    0 1 0;
    sin_AOA 0 cos_AOA];

AOS=user_model.SIDESLIP_ANGLE/180*pi;
cos_AOS=cos(AOS);
sin_AOS=sin(AOS);
rotation_SIDESLIP_ANGLE=[
    cos_AOS -sin_AOS 0;
    sin_AOS cos_AOS 0;
    0 0 1];

free_flow_vector=rotation_AOA*rotation_SIDESLIP_ANGLE*free_flow_vector;
user_model.free_flow_vector=free_flow_vector;

% clear last time record
for monitor_index=1:length(MARKER_MONITERING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITERING(monitor_index),marker_list);
    for element_index=1:marker_list(marker_index).element_number
        marker_element(element_index).streamline_index=[];
    end
end

% for all edge in marker, disturbute opposite edge
for vertex_index=1:length(edge_list)
    edge=edge_list(vertex_index);
    if edge ~= edge_empty
        % means this vertex is not empty
        vertex_ref_list=edge.vertex_ref_list;
        edge_number=edge.edge_number;
        edge.edge_oppo_list=zeros(edge_number,1,'int32');
        
        % add index of edge in HATSEdge
        for edge_index=1:edge.edge_number
            edge.edge_oppo_list(edge_index)=...
                edge_list(vertex_ref_list(edge_index)).getRefIndex(vertex_index);
        end
    end
end

end