classdef HATSElement < handle
    % HATS Element Node Class
    properties
        element_type=[];
        point_index_list=[]; % rank vector
        
        normal_vector=[]; % rank vector
        area=[];
        
        nearby_index_list=[]; % rank vector, contact element list
        nearby_number=int8(0);
    end
    methods
        function HATS_element=HATSElement...
                (element_type,point_index_list)
            HATS_element.element_type=element_type;
            HATS_element.point_index_list=point_index_list;
        end
    end
end