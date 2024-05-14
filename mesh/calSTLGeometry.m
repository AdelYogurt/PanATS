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

% write each marker to file
for marker_index=1:length(marker_name_list)
    marker_name=marker_name_list{marker_index};
    marker=mesh_data.(marker_name);

    element_list=marker.element_list;
    TPL=reshape(element_list',9,[])';

    area=area+sum(vecnorm(cross(TPL(:,4:6)-TPL(:,1:3),TPL(:,7:9)-TPL(:,4:6),2),2,2)/2);
    volume=volume+sum((-TPL(:,7).*TPL(:,5).*TPL(:,3)+TPL(:,4).*TPL(:,8).*TPL(:,3)+TPL(:,7).*TPL(:,2).*TPL(:,6)+...
    -TPL(:,1).*TPL(:,8).*TPL(:,6)-TPL(:,4).*TPL(:,2).*TPL(:,9)+TPL(:,1).*TPL(:,5).*TPL(:,9)))/6;

end

end
