function [area,volume]=solveGeometry()
% calclulate element volume and surface area
%
global user_model

marker_list=user_model.marker_list;
point_list=user_model.point_list;

area=0;
volume=0;

base_z=user_model.min_bou(3)-1;

% calculate each moniter marker
for moniter_index=1:length(user_model.MARKER_MONITORING)
    [marker_element,marker_index]=getMarkerElement...
        (user_model.MARKER_MONITORING(moniter_index),marker_list);
    for element_index=1:marker_list{marker_index,2}
        % area is sum of all element area
        area=area+marker_element(element_index).area;

        % volume is up volume subtract down volume
        pnt_1=point_list(marker_element(element_index).point_index_list(1),1:3);
        pnt_2=point_list(marker_element(element_index).point_index_list(2),1:3);
        pnt_3=point_list(marker_element(element_index).point_index_list(3),1:3);

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

        volume=volume-VL;
    end
end