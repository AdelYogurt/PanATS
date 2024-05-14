function displayMesh(mesh_data,marker_name_list,axe_hdl)
% draw mesh in mesh_data
%
if nargin < 3
    axe_hdl=[];
    if nargin < 2
        marker_name_list=[];
    end
end

if isempty(axe_hdl), axe_hdl=axes(gcf());end

if isempty(marker_name_list)
    marker_name_list=fieldnames(mesh_data);
end
if ischar(marker_name_list)
    marker_name_list={marker_name_list};
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

if isfield(mesh_data,'geometry') && isfield(mesh_data.geometry,'point_list')
    point_list=mesh_data.geometry.point_list;
    dimension=size(point_list,2);
else
    point_list=[];
    dimension=3;
end

for marker_index=1:length(marker_name_list)
    marker_name=marker_name_list{marker_index};
    if strcmp(marker_name,'point_list') || strcmp(marker_name,'name')
        continue;
    end

    % check if part.name exist in part_name_list
    if ~isfield(mesh_data,marker_name)
        continue;
    end

    marker=mesh_data.(marker_name);

    if strcmp(marker.type,'scatter')
        % scatter format data
        hold on;
        if isfield(marker,'Z')
            scatter3(axe_hdl,marker.X,marker.Y,marker.Z);
        else
            scatter(axe_hdl,marker.X,marker.Y);
            dimension=2;
        end
        hold off;
    elseif strcmp(marker.type,'wgs')
        % LaWGS format data
        hold on;
        if isfield(marker,'Z')
            mesh(axe_hdl,marker.X,marker.Y,marker.Z,'FaceColor','none');
        else
            mesh(axe_hdl,marker.X,marker.Y,'FaceColor','none');
            dimension=2;
        end
        hold off;
    elseif strcmp(marker.type,'stl')
        % stl format data
        element_list=marker.element_list;

        % convert stl mesh data structure to matlab mesh data struct
        element_index=reshape(1:size(element_list,1),3,[])';
        hold on;
        trimesh(element_index,element_list(:,1),element_list(:,2),element_list(:,3))
        hold off;
    elseif strcmp(marker.type,'MIXED2') || strcmp(marker.type,'MIXED3')
        data_list=marker.element_list;
        number_list=marker.number_list;
        element_number=length(number_list);

        % convert CGNS mesh data structure to matlab mesh data struct
        element_list=zeros(element_number,3);
        data_idx=0;
        for element_index=1:element_number
            node_num=number_list(element_index);
            if size(element_index,2) < node_num
                element_list=[element_index,nan(elem_num,node_num-size(element_index,2))];
            end
            element_list(element_index,:)=[data_list(data_idx+(1:node_num)+1);data_list(data_idx+2)];

            data_idx=data_idx+node_num+1;
        end
        patch(axe_hdl,'Faces',element_list,'Vertices',point_list,'FaceColor','none')
    else
        element_list=marker.element_list;
        patch(axe_hdl,'Faces',element_list,'Vertices',point_list,'FaceColor','none')
    end
end

axis equal;
xlabel('x');
ylabel('y');

x_range=xlim();
y_range=ylim();
if dimension ~= 2
    zlabel('z');
    view(3);
    z_range=zlim();
    center=[mean(x_range),mean(y_range),mean(z_range)];
    range=max([x_range(2)-x_range(1),y_range(2)-y_range(1),z_range(2)-z_range(1)])/2;
    xlim([center(1)-range,center(1)+range]);
    ylim([center(2)-range,center(2)+range]);
    zlim([center(3)-range,center(3)+range]);
else
    center=[mean(x_range),mean(y_range)];
    range=max([x_range(2)-x_range(1),y_range(2)-y_range(1)])/2;
    xlim([center(1)-range,center(1)+range]);
    ylim([center(2)-range,center(2)+range]);
end

end


