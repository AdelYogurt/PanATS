function preModelStreamline()
% base on half edge data structure to distribute nearby element of element
% method is basing on vertex.neaby_element_list disturbuting nearby element
%
% copyright Adel 2023.03
%
global user_model

vertex_empty=user_model.vertex_empty;
vertex_list=user_model.vertex_list;

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

% initialize set nearby element number to zeros to avoid repeat run error
for moniter_index=1:length(MARKER_MONITERING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITERING{moniter_index},marker_list);
    for element_index=1:marker_list(marker_index).element_number
        marker_element(element_index).element_nearby_number=0;
    end
end

% for all vertex in marker, disturbute nearby element
for vertex_index=1:length(vertex_list)
    vertex=vertex_list(vertex_index);
    if vertex ~= vertex_empty
        % means this vertex is not empty
        % add reference element into dual edge element
        for nearby_index=1:vertex.nearby_number
            point_nearby_index=vertex.point_next_list(nearby_index);

            % cheak dual element by hash table(vertex_list)
            vertex_next=vertex_list(point_nearby_index);
            if vertex ~= vertex_empty
                % search vertex_index place in vertex_next.point_next_list
                [exist_flag,ref_index]=judgeMatExistNum...
                    (vertex_next.point_next_list,vertex_index);
                if exist_flag
                    element_dual=vertex_next.element_ref_list(ref_index);
                    element_dual.element_nearby_number=element_dual.element_nearby_number+1;
                    element_dual.element_nearby_list(element_dual.element_nearby_number)=...
                        vertex.element_ref_list(nearby_index);
                end
            end
        end
    end
end

% for all moniter marker element, generate element_outflow_boolean
for moniter_index=1:length(MARKER_MONITERING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITERING{moniter_index},marker_list);
    for element_index=1:marker_list(marker_index).element_number
        element=marker_element(element_index);
        element_nearby_number=element.element_nearby_number;
        surface_flow=element.surface_flow;
        element.element_outflow_boolean=true(element_nearby_number,1); % default is outflow element

        for nearby_index=1:element_nearby_number
            % base on 
            if dot(element.center_point-element.element_nearby_list(nearby_index).center_point,...
                    free_flow_vector) > 0
                element.element_outflow_boolean(nearby_index)=true(1);
            end
        end
    end
end

end