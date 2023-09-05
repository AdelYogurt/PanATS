function displayMarker(marker_name_list)
% draw marker point
%
global user_model
if nargin < 1
    marker_name_list=[];
end

geometry=user_model.geometry;
marker_list=user_model.marker_list;
element_list=user_model.element_list;

dimension=geometry.dimension;
point_list=geometry.point_list;

if isempty(marker_name_list)
    marker_name_list=fieldnames(marker_list);
end

if ischar(marker_name_list)
    marker_name_list={marker_name_list};
end

hold on;
for marker_index=1:length(marker_name_list)
    marker_name=marker_name_list{marker_index};
    Elem_idx=marker_list.(marker_name);

    tri_elem_list=zeros(length(Elem_idx),3);
    for elem_idx=1:length(Elem_idx)
        elem=element_list(Elem_idx(elem_idx));
        % element point index
        tri_elem_list(elem_idx,:)=elem.Node_idx(1:3);
    end
    trimesh(tri_elem_list,point_list(:,1),point_list(:,2),point_list(:,3),'FaceColor','None')
end
hold off;

x_range=xlim();
y_range=ylim();
z_range=zlim();

center=[mean(x_range),mean(y_range),mean(z_range)];
range=max([x_range(2)-x_range(1),y_range(2)-y_range(1),z_range(2)-z_range(1)])/2;

axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
view(3);

xlim([center(1)-range,center(1)+range]);
ylim([center(2)-range,center(2)+range]);
zlim([center(3)-range,center(3)+range]);
drawnow;
end
