function displayMarker(marker_name_list)
% draw marker point
%
global user_model
if nargin < 1
    marker_name_list=[];
end

dimension=user_model.geometry.dimension;
marker_list=user_model.marker_list;
point_list=user_model.point_list;

if isempty(marker_name_list)
    marker_name_list={marker_list.name};
end

if ischar(marker_name_list)
    marker_name_list={marker_name_list};
end

for marker_index=1:length(marker_name_list)
    marker_element=getMarkerElement(marker_name_list{marker_index},marker_list);
    for element_index=1:size(marker_element,1)
        element=marker_element(element_index);
        % element point index
        node_index=[element.point_index_list,...
            element.point_index_list(1)];
        line(point_list(node_index,1),point_list(node_index,2),point_list(node_index,3));
        if ~isempty(element.attachment)
            if element.attachment
                patch(point_list(node_index,1),point_list(node_index,2),point_list(node_index,3),'g');
            end

        end
    end
end

axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view(3);

end