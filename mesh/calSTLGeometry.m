function [area,volume]=calSTLGeometry(part)
% calclulate stl mesh volume and surface area
%
mesh_list=part.mesh_list;
mesh_number=length(mesh_list);

area=0;
volume=0;

base_z=0;

for mesh_index=1:mesh_number
    % calculate each mesh
    mesh=mesh_list{mesh_index};

    mesh_element_number=mesh.element_number;
    mesh_element_list=mesh.element_list;

    for element_index=1:mesh_element_number
        % volume is up volume subtract down volume
        pnt_1=mesh_element_list(3*element_index-2,1:3);
        pnt_2=mesh_element_list(3*element_index-1,1:3);
        pnt_3=mesh_element_list(3*element_index,1:3);

        % area is sum of all element area
        area=area+norm(cross(pnt_2-pnt_1,pnt_3-pnt_2))/2;

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
end