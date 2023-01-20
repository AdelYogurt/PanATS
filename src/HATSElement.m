classdef HATSElement < handle
    % HATS Element Node Class
    properties
        % base properties
        element_type=[];
        point_index_list=[]; % rank vector

        marker_index=[];
        
        element_nearby_list=[]; % contact element list, HATSElement pointer
        element_nearby_number=int8(0);

        % geometry properties
        normal_vector=[];
        area=[];
        center_point=[];
    end
    methods
        function HATS_element=HATSElement...
                (element_type,point_index_list)
            HATS_element.element_type=element_type;
            HATS_element.point_index_list=point_index_list;
        end
    end
end