function solveModelStreamline()
% calculate element reference streaamline length
%
% reference:
% [1]KENWRIGHT D N. Automatic detection of open and closed
% separation and attachment lines; proceedings of the Proceedings
% Visualization '98 (Cat No98CB36276), F 18-23 Oct. 1998, 1998 [C].
%
global g_model

dimension=g_model.dimension;
point_list=g_model.point_list;
marker_list=g_model.marker_list;
MARKER_MONITERING=g_model.MARKER_MONITORING;

geometry_torlance=1e-12;

% calculate inflow vector
free_flow_vector=[1;0;0];

AOA=g_model.AOA/180*pi;
cos_AOA=cos(AOA);
sin_AOA=sin(AOA);
rotation_AOA=[
    cos_AOA 0 -sin_AOA;
    0 1 0;
    sin_AOA 0 cos_AOA];

AOS=g_model.SIDESLIP_ANGLE/180*pi;
cos_AOS=cos(AOS);
sin_AOS=sin(AOS);
rotation_SIDESLIP_ANGLE=[
    cos_AOS -sin_AOS 0;
    sin_AOS cos_AOS 0;
    0 0 1];

free_flow_vector=rotation_AOA*rotation_SIDESLIP_ANGLE*free_flow_vector;
g_model.free_flow_vector=free_flow_vector;

% calculate surface velocity
surface_flow_list=[zeros(size(point_list,1),4)]; % data, repeat times
for moniter_index=1:length(MARKER_MONITERING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITERING{moniter_index},marker_list);
    for element_index=1:marker_list(marker_index).element_number
        element=marker_element(element_index);
        normal_vector=element.normal_vector;
        point_index_list=element.point_index_list;
        dot_nor_vec=dot(free_flow_vector,normal_vector);
        Ve=(free_flow_vector'-normal_vector*dot_nor_vec);
        element.surface_flow=Ve;

        % add data to each point
        surface_flow_list(point_index_list,1:3)=...
            surface_flow_list(point_index_list,1:3)+Ve;
        surface_flow_list(point_index_list,4)=...
            surface_flow_list(point_index_list,4)+1;
        
        % process symmetry
        if element.element_nearby_number < length(element.point_index_list)
            % means element on symmetry face
            switch g_model.SYMMETRY
                case 'XOY'
                    flip_operate=[1,1,-1];
                    identify_dimension=3;
                case 'YOZ'
                    flip_operate=[-1,1,1];
                    identify_dimension=1;
                case 'ZOX'
                    flip_operate=[1,-1,1];
                    identify_dimension=2;
                otherwise
                    error('solveModelStreamline: nuknown SYMMETRY type');
            end
            element_symmetry=HATSElement(element.element_type,...
                element.point_index_list);
            element_symmetry.center_point=element.center_point.*...
                flip_operate;
            element_symmetry.surface_flow=element.surface_flow.*...
                flip_operate;
            element_symmetry.normal_vector=element.normal_vector.*...
                flip_operate;

            % add to element_nearby_list
            element.element_nearby_number=element.element_nearby_number+1;
            element.element_nearby_list(element.element_nearby_number)=...
                element_symmetry;

            % identify which point of element on symmetry face
            % add symmetry velocity to point
            Ve(identify_dimension)=-Ve(identify_dimension);
            for node_index=1:lenght(point_index_list)
                point_index=point_index_list(node_index);
                if abs(point_list(point_index,dimension)) <= ...
                        geometry_torlance
                    % on symmetry face
                    % add data to each point
                    surface_flow_list(point_index,1:3)=...
                        surface_flow_list(point_index,1:3)+Ve;
                    surface_flow_list(point_index,4)=...
                        surface_flow_list(point_index,4)+1;
                end
            end
        end
    end
end
% average velocity
surface_flow_list(point_index_list,1:3)=...
    surface_flow_list(point_index_list,1:3)./surface_flow_list(point_index_list,4);

% identify stagnation element
% base on reference [1] to extract
for moniter_index=1:length(MARKER_MONITERING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITERING{moniter_index},marker_list);
    for element_index=1:marker_list(marker_index).element_number
        % base on around element to determine if stagnation element
        element=marker_element(element_index);
        if element.normal_vector*free_flow_vector > 0
            element.stagnation=int8(0);
            continue;
        end
        abs_Ve=norm(element.surface_flow);
        ez=element.normal_vector;
        if abs(abs_Ve) == 0
            element.stagnation=int8(0);
            continue;
        end

        % only consider upwind element have surface flow
        ex=element.surface_flow/abs_Ve;
        ey=cross(ez,ex);

        element_nearby_number=element.element_nearby_number;
        element_nearby_list=element.element_nearby_list(1:element_nearby_number);
        element_point_list=point_list(element.point_index_list,:);

        if element.normal_vector*free_flow_vector < -0.98
            surface_flow=element.surface_flow;
            quiver3(0,0,0,surface_flow(1),surface_flow(2),surface_flow(3),0.001);
            drawElement(element,point_list,element.center_point)
            hold on;
            for nearby_index=1:element_nearby_number
                drawFlow(element_nearby_list(nearby_index),...
                    element_nearby_list(nearby_index).center_point-element.center_point)
                drawElement(element_nearby_list(nearby_index),point_list,element.center_point)
            end
            hold off;
            view(3);
            xlabel('x');
            ylabel('y');
            zlabel('z');
            axis equal;
            disp('head');
        end

        % all vector just simply project to element panel
        Ve_list=zeros(element_nearby_number,3);
        d_center_point_list=zeros(element_nearby_number,3);
        for nearby_index=1:element_nearby_number
            Ve_list(nearby_index,:)=[element_nearby_list(nearby_index).surface_flow]*[ex;ey;ez]';
            d_center_point_list(nearby_index,:)=(element_nearby_list(nearby_index).center_point-...
                element.center_point)*[ex;ey;ez]';
        end
        d_edge_list=(element_point_list-element.center_point)*[ex;ey;ez]';
        d_center_point_list(:,end)=1;

        % fit flow
        matrix=reshape([d_center_point_list,zeros(element_nearby_number,6),d_center_point_list]',6,element_nearby_number*2)';
        matrix=[0,0,1,0,0,0;0,0,0,1,0,0;matrix];
        Ve_list=[abs_Ve,0;Ve_list(:,1:2)]';
        colume=Ve_list(:);

        coefficient=matrix\colume;

        A=[coefficient(3);coefficient(6)];
        BC=[coefficient(1:2)';coefficient(4:5)'];

        if rank(BC) == 2
            eig_value=eig(BC);
            stagnation_point=-BC\A;

            if sum(eig_value > 0)
                % judge if stagnation_point is in element
                d_edge_list=d_edge_list(:,1:2)-stagnation_point';
                cross_time=0;
                for edge_index=1:size(d_edge_list,1)
                    edge_index_next=edge_index+1;
                    if edge_index_next > size(d_edge_list,1)
                        edge_index_next=1;
                    end

                    if ((-d_edge_list(edge_index_next,2)) * ...
                            ((d_edge_list(edge_index,2)))) < -eps
                        % if cross point not in line y range
                        continue;
                    end

                    A=d_edge_list(edge_index,2)-d_edge_list(edge_index_next,2);
                    C=d_edge_list(edge_index_next,2)*d_edge_list(edge_index,1)-...
                        d_edge_list(edge_index_next,1)*d_edge_list(edge_index,2);
                    cross_point=-C/A;
                    if abs(A) > 0
                        if cross_point > 0
                            if ((cross_point - d_edge_list(edge_index_next,1)) * ...
                                    ((d_edge_list(edge_index,1)) - cross_point)) >= -eps
                                % cross point must in line range
                                cross_time=cross_time+1;
                            end
                        end
                    else
                        if ((-d_edge_list(edge_index_next,1)) * ...
                                ((d_edge_list(edge_index,1)))) >=0
                            % cross point must in line range
                            cross_time=cross_time+1;
                        end
                    end
                end
                if cross_time == 1
                    element.stagnation=int8(1);
                else
                    element.stagnation=int8(0);
                end
            else
                element.stagnation=int8(0);
            end
        else
            element.stagnation=int8(0);
        end
    end
end

% % calculate streamline length
% inflow_direction_list=zeros(marker_element_number,3);
% streamline_length_list=zeros(marker_element_number,1);
%
% bottom_index_list=[];
% for element_index=1:marker_element_number
%     center_point=ADtree_marker_element.center_point_list(element_index,:);
%     out_index_list=HATS_element_list(element_index).out_index_list;
%     normal_vector=HATS_element_list(element_index).normal_vector;
%
%     if abs(normal_vector*free_flow_vector-1) < geometry_torlance
%         bottom_index_list=[bottom_index_list,element_index];
%         continue;
%     end
%
%     inflow_direction=inflow_direction_list(element_index,:);
%     streamline_length=streamline_length_list(element_index);
%
%     for out_node_index__=1:length(out_index_list)
%         out_node_index=out_index_list(out_node_index__);
%
%         out_center_point=ADtree_marker_element.center_point_list(out_node_index,:);
%
%         % base on normal_vector calculate inflow_direction and streamline_length
%         inflow_direction_out=(out_center_point-center_point);
%         streamline_length_out=sqrt(sum(inflow_direction_out.^2));
%         inflow_direction_out=inflow_direction_out/streamline_length_out;
%
%         % correct length and vector base on parent inflow_direction and streamline_length
%         if ~isempty(inflow_direction)
%             inflow_direction_out=...
%                 inflow_direction_out*streamline_length_out+...
%                 inflow_direction*streamline_length; % vector plus
%             streamline_length_out=sqrt(sum(inflow_direction_out.^2));
%             inflow_direction_out=inflow_direction_out/streamline_length_out;
%         end
%
%         % compare exist steamline, if shorter than repace
%         if (streamline_length_list(out_node_index) ~= 0)
%             if streamline_length_list(out_node_index) > streamline_length_out
%                 inflow_direction_list(out_node_index,:)=inflow_direction_out;
%                 streamline_length_list(out_node_index)=streamline_length_out;
%             end
%         else
%             inflow_direction_list(out_node_index,:)=inflow_direction_out;
%             streamline_length_list(out_node_index)=streamline_length_out;
%         end
%     end
% end
% inflow_direction_list(bottom_index_list,:)=repmat(free_flow_vector',length(bottom_index_list),1);
% streamline_length_list(bottom_index_list)=max(streamline_length_list);
%
% streamline_output.inflow_direction_list=inflow_direction_list;
% streamline_output.streamline_length_list=streamline_length_list;
end
function drawElement(element,point_list,origin_point)
element_point_list=point_list(element.point_index_list,:);
element_point_list=(element_point_list-origin_point);
order=[1,2,3,1];
line(element_point_list(order,1),element_point_list(order,2),element_point_list(order,3));
end
function drawFlow(element,point)
surface_flow=element.surface_flow;
quiver3(point(1),point(2),point(3),surface_flow(1),surface_flow(2),surface_flow(3),0.001);
end