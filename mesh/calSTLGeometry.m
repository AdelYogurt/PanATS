function [area,volume]=calSTLGeometry(mesh_data,marker_name_list)
% calclulate stl mesh volume and surface area
%
% input:
% mesh_data, marker_name_list(default all markers)
%
% output:
% area, volume
%
if nargin < 2
    marker_name_list=[];
end
if isempty(marker_name_list)
    marker_name_list=fieldnames(mesh_data);
end
marker_index=1;
while marker_index <= length(marker_name_list)
    marker_name=marker_name_list{marker_index};
    if strcmp(marker_name,'geometry')
        marker_name_list(marker_index)=[];
    else
        marker_index=marker_index+1;
    end
end

area=0;
volume=0;
base_z=0;

% write each marker to file
for marker_index=1:length(marker_name_list)
    marker_name=marker_name_list{marker_index};
    marker=mesh_data.(marker_name);

    element_list=marker.element_list;
    element_number=size(element_list,1)/3;

    for element_index=1:element_number
        % volume is up volume subtract down volume
        pnt_1=element_list(3*element_index-2,1:3);
        pnt_2=element_list(3*element_index-1,1:3);
        pnt_3=element_list(3*element_index,1:3);

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