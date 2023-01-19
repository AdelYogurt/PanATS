function displayMarker(marker_name_list)
% draw marker point
%
global g_model
if nargin < 1
    marker_name_list=[];
end

dimension=g_model.dimension;
marker_list=g_model.marker_list;
point_list=g_model.point_list;

if isempty(marker_name_list)
    marker_name_list=marker_list(:,1);
else
    if ischar(marker_name_list)
        marker_name_list={marker_name_list};
    end
end

if dimension == 2
    %     line(point_list(point_index_list+1,1),point_list(point_index_list+1,2),...
    %         'LineStyle','none','Marker','o');
    for marker_index=1:length(marker_name_list)
        marker_element=getMarkerElement(marker_name_list{marker_index});
        for element_index=1:size(marker_element,1)
            element=marker_element(element_index,:);
            element_type=element(1);

            % element point index
            switch(element_type)
                case 3
                    point_number=2;
                case 5
                    point_number=3;
                case 9
                    point_number=4;
                case 10
                    point_number=4;
                case 12
                    point_number=8;
                case 13
                    point_number=6;
                case 14
                    point_number=5;
            end

            node_index=[element(2:1+point_number),element(2)];
            line(point_list(node_index,1),point_list(node_index,2));
        end
    end
elseif dimension == 3
    %     line(point_list(point_index_list+1,1),point_list(point_index_list+1,2),point_list(point_index_list+1,3),...
    %         'LineStyle','none','Marker','o');
    for marker_index=1:length(marker_name_list)
        marker_element=getMarkerElement(marker_name_list{marker_index});
        for element_index=1:size(marker_element,1)
            element=marker_element(element_index,:);
            element_type=element(1);

            % element point index
            switch(element_type)
                case 3
                    point_number=2;
                case 5
                    point_number=3;
                case 9
                    point_number=4;
                case 10
                    point_number=4;
                case 12
                    point_number=8;
                case 13
                    point_number=6;
                case 14
                    point_number=5;
            end

            node_index=[element(2:1+point_number),element(2)];
            line(point_list(node_index,1),point_list(node_index,2),point_list(node_index,3));
        end
    end

    zlabel('z');
    view(3);
end
axis equal;
xlabel('x');
ylabel('y');

    function Marker_Element=getMarkerElement(marker_name)
        % return specified marker
        % marker is include all element unit of marker
        % element include element type and point index
        %
        Marker_Element=[];
        for marker_index=1:size(marker_list,1)
            if strcmp(marker_list{marker_index,1},marker_name)
                Marker_Element=marker_list{marker_index,3};
            end
        end
        if isempty(Marker_Element)
            error('getMarkerElement: no marker found')
        end
    end
end