classdef HATSVertex < handle
    % half edge structure
    % notice point(self) to nearby point and element(nearby) comply with the right hand rule
    properties
        element_ref_list=[]; % pointer to nearby element form a halfedge, this sort face of halfedge
        point_next_list=[]; % start from 4, if no enough will increase
        point_prev_list=[]; % start from 4, if no enough will increase
        nearby_number=[];
    end
end