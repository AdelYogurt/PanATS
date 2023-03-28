classdef HATSElement < handle
    % HATS Element Node Class
    properties
        % base properties
        element_type=[];
        point_index_list=[]; % rank vector

        marker_index=[];
        
        element_nearby_list=[]; % contact element list, HATSElement pointer
        element_nearby_number=int8(0);
        element_outflow_boolean=[]; % if nearby element is outflow of self element

        % geometry properties
        normal_vector=[];
        area=[];
        center_point=[];
        surface_flow=[]; % norm range is -1 to 1, 1 means equal to V_1
        streamline_length=[];
        inflow_vector=[];
        attachment=[];
    end
    methods
        function HATS_element=HATSElement...
                (element_type,point_index_list)
            HATS_element.element_type=element_type;
            HATS_element.point_index_list=point_index_list;
        end
    end
end