function [area,area_x,area_y,area_z,volume]=solveGeometry()
% calclulate element volume and surface area
%
% copyright Adel 2023.03
%
global user_model

config=user_model.config;
geometry=user_model.geometry;
element_list=user_model.element_list;
SYMMETRY=config.SYMMETRY;

% load geometry
dimension=geometry.dimension;
point_list=geometry.point_list;
area_list=geometry.area_list;
normal_vector_list=geometry.normal_vector_list;

% area is sum of all element area
area=sum(area_list);
area_x=sum(area_list.*max(normal_vector_list(:,1),0));
area_y=sum(area_list.*max(normal_vector_list(:,2),0));
area_z=sum(area_list.*max(normal_vector_list(:,3),0));

% load all point list
elem_num=length(element_list);
TPL=zeros(elem_num,9);
tri_num=0;
for elem_idx=1:elem_num
    Node_idx=element_list(elem_idx).Node_idx;
    node_num=element_list(elem_idx).node_num;

    % volume is up volume subtract down volume
    for node_idx=2:node_num-1
        pnt_1=point_list(Node_idx(1),1:3);
        pnt_2=point_list(Node_idx(node_idx),1:3);
        pnt_3=point_list(Node_idx(node_idx+1),1:3);
        tri_num=tri_num+1;

        if size(TPL,1) < tri_num
            TPL=[TPL;zeros(elem_num,9)];
        else
            TPL(tri_num,:)=[pnt_1,pnt_2,pnt_3];
        end
    end
end
TPL(tri_num,:)=[];

volume=sum((-TPL(:,7).*TPL(:,5).*TPL(:,3)+TPL(:,4).*TPL(:,8).*TPL(:,3)+TPL(:,7).*TPL(:,2).*TPL(:,6)+...
    -TPL(:,1).*TPL(:,8).*TPL(:,6)-TPL(:,4).*TPL(:,2).*TPL(:,9)+TPL(:,1).*TPL(:,5).*TPL(:,9)))/6;

if SYMMETRY
    area=area*2;
    volume=volume*2;
    switch SYMMETRY
        case 'YOZ'
            area_y=area_y*2;
            area_z=area_z*2;
        case 'ZOX'
            area_z=area_z*2;
            area_x=area_x*2;
        case 'XOY'
            area_x=area_x*2;
            area_y=area_y*2;
        otherwise
            error('solveGeometry: nuknown SYMMETRY type');
    end
end

    function VL=calProjVolume(pnt_1,pnt_2,pnt_3)
        % calculate volume base on three point
        %

        VF=det(...
            [1 pnt_1(1) pnt_1(2) pnt_1(3);
            1 pnt_2(1) pnt_2(2) pnt_2(3);
            1 pnt_3(1) pnt_3(2) pnt_3(3);
            1 pnt_1(1) pnt_1(2) 0;]);
        VS=det(...
            [1 pnt_2(1) pnt_2(2) pnt_2(3);
            1 pnt_2(1) pnt_2(2) 0;
            1 pnt_3(1) pnt_3(2) pnt_3(3);
            1 pnt_1(1) pnt_1(2) 0;]);
        VT=det(...
            [1 pnt_1(1) pnt_1(2) 0;
            1 pnt_3(1) pnt_3(2) 0;
            1 pnt_2(1) pnt_2(2) 0;
            1 pnt_3(1) pnt_3(2) pnt_3(3);]);
        VL=-(VF+VS+VT)/6;
    end
end
