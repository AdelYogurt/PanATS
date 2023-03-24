function preModelStreamline()
% base on half edge data structure to distribute nearby element of element
% method is basing on vertex.neaby_element_list disturbuting nearby element
%
% copyright Adel 2023.03
%
global user_model

vertex_empty=user_model.vertex_empty;
element_empty=user_model.element_empty;
vertex_list=user_model.vertex_list;

marker_list=user_model.marker_list;
MARKER_MONITERING=user_model.MARKER_MONITORING;

% initialize set nearby element number to zeros to avoid repeat run error
for moniter_index=1:length(MARKER_MONITERING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITERING{moniter_index},marker_list);
    for element_index=1:marker_list(marker_index).element_number
        marker_element(element_index).element_nearby_number=0;
    end
end

% for all vertex in marker, disturbute nearby element
for point_base_index=1:length(vertex_list)
    vertex=vertex_list(point_base_index);
    if vertex ~= vertex_empty
        % means this vertex is not empty
        % add this element into dual edge element
        for nearby_index=1:vertex.nearby_number
            point_nearby_index=vertex.point_next_list(nearby_index);
            % cheak dual element by hash table(vertex_list)
            vertex_dual=vertex_list(point_nearby_index);
            if vertex ~= vertex_empty
                [exist_flag,index]=judgeMatExistNum...
                    (vertex_dual.point_next_list,point_base_index);
                if exist_flag
                    element_dual=vertex_dual.element_ref_list(index);
                    element_dual.element_nearby_number=element_dual.element_nearby_number+1;
                    element_dual.element_nearby_list(element_dual.element_nearby_number)=...
                        vertex.element_ref_list(nearby_index);
                end
            end
        end
    end
end

end