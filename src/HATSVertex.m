classdef HATSVertex < handle
    properties
        point_nearby_list=[]; % start from 4, if no enough will increase
        element_list=[]; % pointer to nearby element form a halfedge, this sort face of halfedge
        nearby_number=[];
    end
end