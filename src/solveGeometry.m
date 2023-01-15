function [area,volume]=solveGeometry()
% calclulate element volume and surface area
%
global g_geometry g_Point g_Element g_Marker...
    ADtree_marker_element HATS_element_list

% area is sum of all element area
area=sum([HATS_element_list.area]);

% volume is up volume subtract down volume
base_z=g_geometry.min_bou(3)-1;
volume=0;
for element_index=1:length(HATS_element_list)
    pnt_1=g_Point(HATS_element_list(element_index).point_index_list(1),1:3);
    pnt_2=g_Point(HATS_element_list(element_index).point_index_list(2),1:3);
    pnt_3=g_Point(HATS_element_list(element_index).point_index_list(3),1:3);
    
    % element volume to base plate
    VF=det(...
        [1 pnt_1(1) pnt_1(2) pnt_1(3);
        1 pnt_2(1) pnt_2(2) pnt_2(3);
        1 pnt_3(1) pnt_3(2) pnt_3(3);
        1 pnt_1(1) pnt_1(2) base_z;]);
    VS=det(...
        [1 pnt_2(1) pnt_2(2) pnt_2(3);
        1 pnt_2(1) pnt_2(2) base_z;
        1 pnt_3(1) pnt_3(2) pnt_3(3);
        1 pnt_1(1) pnt_1(2) base_z;]);
    VT=det(...
        [1 pnt_1(1) pnt_1(2) base_z;
        1 pnt_3(1) pnt_3(2) base_z;
        1 pnt_2(1) pnt_2(2) base_z;
        1 pnt_3(1) pnt_3(2) pnt_3(3);]);
    
    VL=(VF+VS+VT)/6;
    
    if HATS_element_list(element_index).normal_vector(3) > 0
        volume=volume-VL;
    else
        volume=volume-VL;
    end
end
end