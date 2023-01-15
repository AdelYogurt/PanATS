classdef HATSElement < handle
    % HATS Element Node Class
    properties
        element_type;
        point_index_list; % rank vector
        
        normal_vector; % rank vector
        area;
        
        in_index_list; % rank vector, inflow element index
        out_index_list; % rank vector, outflow element index
    end
    methods
        function HATS_element=HATSElement...
                (element_type,point_index_list)
            HATS_element.element_type=element_type;
            HATS_element.point_index_list=point_index_list;
        end
    end
end