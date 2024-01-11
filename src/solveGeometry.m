function [area,volume]=solveGeometry()
% calclulate element volume and surface area
%
% copyright Adel 2023.03
%
global user_model

geometry=user_model.geometry;

dimension=geometry.dimension;
point_list=geometry.point_list;
element_list=user_model.element_list;

SYMMETRY=user_model.SYMMETRY;

% load geometry
area_list=geometry.area_list;

area=0;
volume=0;
base_z=0;

elem_num=length(element_list);
for elem_idx=1:elem_num
    Node_idx=element_list(elem_idx).Node_idx;
    node_num=element_list.node_num;

    % area is sum of all element area
    area=area+area_list(elem_idx);

    % volume is up volume subtract down volume
    for node_idx=2:node_num-1
        pnt_1=point_list(Node_idx(1),1:3);
        pnt_2=point_list(Node_idx(node_idx),1:3);
        pnt_3=point_list(Node_idx(node_idx+1),1:3);
        VL=calProjVolume(pnt_1,pnt_2,pnt_3,base_z);
    end

    volume=volume+VL;
end

if SYMMETRY
    area=area*2;
    volume=volume*2;
end

    function VL=calProjVolume(pnt_1,pnt_2,pnt_3,base_z)
        % calculate volume base on three point
        %
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
        VL=-(VF+VS+VT)/6;
    end
end
