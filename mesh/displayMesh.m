function displayMesh(grid,marker_name_list)
% draw mesh in grid
%
if nargin < 2
    marker_name_list=[];
end

% default draw all marker
if isempty(marker_name_list)
    marker_name_list=fieldnames(grid);
end

if ischar(marker_name_list)
    marker_name_list={marker_name_list};
end

if isfield(grid,'point_list')
    point_list=grid.point_list;
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
    if ~isfield(grid,marker_name)
        continue;
    end

    marker=grid.(marker_name);

    if strcmp(marker.type,'wgs')
        % LaWGS format data
        hold on;
        surface(marker.X,marker.Y,marker.Z,'FaceColor','none');
        hold off;
    elseif strcmp(marker.type,'stl')
        % stl format data
        element_list=marker.element_list;
        element_number=size(element_list,1)/3;

        for element_index=1:element_number
            point_index_list=[3*element_index-2:3*element_index,3*element_index-2];
            line(element_list(point_index_list,1), ...
                element_list(point_index_list,2), ...
                element_list(point_index_list,3))
        end
    elseif strcmp(marker.type,'MIXED2') || strcmp(marker.type,'MIXED3')
        element_list=marker.element_list;
        number_list=marker.number_list;
        element_number=length(number_list);

        node_index=0;
        for element_index=1:element_number
            node_number=number_list(element_index);

            % element point index
            point_index_list=[element_list(node_index+(1:node_number));element_list(node_index+1)];

            if dimension == 2
                line(point_list(point_index_list,1), ...
                    point_list(point_index_list,2));
            else
                line(point_list(point_index_list,1), ...
                    point_list(point_index_list,2), ...
                    point_list(point_index_list,3));
            end

            node_index=node_index+node_number;
        end
    else
        element_list=marker.element_list;
        element_number=size(element_list,1);

        for element_index=1:element_number
            % element point index
            point_index_list=[element_list(element_index,:),element_list(element_index,1)];

            if dimension == 2
                line(point_list(point_index_list,1), ...
                    point_list(point_index_list,2));
            else
                line(point_list(point_index_list,1), ...
                    point_list(point_index_list,2), ...
                    point_list(point_index_list,3));
            end
        end
    end
end

axis equal;
xlabel('x');
ylabel('y');

if dimension ~= 2
    zlabel('z');
    view(3);
end

end


