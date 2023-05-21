function displayPart(part_list,point_list,part_name_list)
% draw part list
%
if nargin < 3
    part_name_list=[];
    if nargin < 2
        point_list=[];
    end
end

if ~iscell(part_list)
    part_list={part_list};
end

part_number=length(part_list);

% default draw all part
if isempty(part_name_list)
    part_name_list=cell(1,part_number);
    for part_index=1:part_number
        part_name_list{part_index}=part_list{part_index}.name;
    end
end

if ischar(part_name_list)
    part_name_list={part_name_list};
end

for part_index=1:part_number
    part=part_list{part_index};

    % check if part.name exist in part_name_list
    exist=false(1);
    for draw_index=1:length(part_name_list)
        if strcmp(part.name,part_name_list{draw_index})
            exist=true(1);
            break;
        end
    end

    % if not exist, continue next
    if ~exist
        continue;
    end

    mesh_list=part.mesh_list;

    for mesh_index=1:length(mesh_list)
        mesh=mesh_list{mesh_index};

        if strcmp(mesh.element_type,'wgs')
            % LaWGS format data
            hold on;
            surf(mesh.X,mesh.Y,mesh.Z,'FaceColor','none');
            hold off;
        elseif strcmp(mesh.element_type,'stl')
            % stl format data
            element_list=mesh.element_list;
            element_number=size(element_list,1)/3;

            for element_index=1:element_number
                point_index_list=[3*element_index-2:3*element_index,3*element_index-2];
                line(element_list(point_index_list,1), ...
                    element_list(point_index_list,2), ...
                    element_list(point_index_list,3))
            end
        else
            element_list=mesh.element_list;
            element_number=size(element_list,1);

            switch mesh.element_ID
                case 5 % triangle element
                    point_number=3;
                case 9 % quadrilateral element
                    point_number=4;
            end

            for element_index=1:element_number
                % element point index
                point_index_list=[element_list(element_index,1:point_number),...
                    element_list(element_index,1)];
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
zlabel('z');
view(3);

end
