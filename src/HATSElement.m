classdef HATSElement < handle
    % HATS Element Node Class
    properties
        % base properties
        element_type=[]; % int8
        % rank vector, also implicit storage edge index
        point_index_list=[]; % int32
        % marker_index
        marker_index=[]; % int8
        % edge index in HATSEdge
        edge_ref_list=[]; % int8
    end
    properties % geometry data
        normal_vector=[];
        area=[];
        center_point=[];
    end
    properties % streamline data
        % norm range is -1 to 1, 1 means equal to V_1
        surface_flow=[];
        % if this element is attachment line element of stagnation element
        attachment=[]; % boolean
        % index of streamline which across this element
        streamline_ref=[]; % streamline reference, index of segment in streamline
    end
    methods
        function HATS_element=HATSElement...
                (element_type,point_index_list)
            HATS_element.element_type=element_type;
            HATS_element.point_index_list=point_index_list;
            HATS_element.edge_ref_list=zeros(length(point_index_list),1,'int8');
        end
    end
end