function displayMesh(axe_hdl,mesh_data,marker_name_list,option)
% draw mesh in mesh_data
%
if nargin < 4
    option=[];
    if nargin < 3
        marker_name_list=[];
    end
end
if isempty(axe_hdl),axe_hdl=gca();end

patch_option=struct('FaceColor','none');
if ~isempty(option)
    names=fieldnames(option);
    for idx=1:length(names),patch_option.(names{idx})=option.(names{idx});end
end

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

hold on;
for marker_index=1:length(marker_name_list)
    marker_name=marker_name_list{marker_index};

    % check if part.name exist in part_name_list
    if ~isfield(mesh_data,marker_name),continue;end

    marker=mesh_data.(marker_name);

    if strcmp(marker.type,'scatter')
        % scatter format data
        if isfield(marker,'Z')
            scatter3(axe_hdl,marker.X,marker.Y,marker.Z,patch_option);
        else
            scatter(axe_hdl,marker.X,marker.Y,patch_option);
            dimension=2;
        end
    elseif strcmp(marker.type,'wgs')
        % LaWGS format data
        if isfield(marker,'Z')
            surface(axe_hdl,marker.X,marker.Y,marker.Z,patch_option);
        else
            surface(axe_hdl,marker.X,marker.Y,patch_option);
            dimension=2;
        end
    elseif strcmp(marker.type,'stl')
        % stl format data
        element_list=marker.element_list;
        point_list_marker=element_list;
        element_list_marker=reshape(1:size(element_list,1),3,[])';
        
        patch(axe_hdl,'Faces',element_list_marker,'Vertices',point_list_marker,patch_option)
    elseif strcmp(marker.type,'MIXED2') || strcmp(marker.type,'MIXED3')
        element_list=marker.element_list;
        number_list=marker.number_list;
        element_number=length(number_list);

        % convert CGNS mix mesh data structure to matlab mesh data struct
        element_list_marker=zeros(element_number,3);node_idx=1;
        for element_index=1:element_number
            node_num=number_list(element_index);
            node_index_list=element_list(node_idx+(1:node_num));
            if size(element_list_marker,2) < node_num
                element_list_marker=[element_list_marker,nan(element_number,node_num-size(element_list_marker,2))];
            end
            element_list_marker(element_index,:)=node_index_list;
            node_idx=node_idx+node_num+1;
        end

        patch(axe_hdl,'Faces',element_list_marker,'Vertices',point_list,patch_option)
    else
        element_list=marker.element_list;
        patch(axe_hdl,'Faces',element_list,'Vertices',point_list,patch_option);
    end
end
hold off;

axis equal;
xlabel('x');
ylabel('y');

x_range=xlim();
y_range=ylim();
if dimension ~= 2
    zlabel('z');view(3);
    
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


