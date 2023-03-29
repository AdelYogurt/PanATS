function preModelPanel()
% function to prepare model for panel method
% including:
% generate edge_list, construct half of edge data structure,
% and make marker monitered element normal vector consisted
% calculate monitor element geomertry, normal vector and area of element
%
% copyright Adel 2023.03
%
global user_model

point_list=user_model.point_list;
marker_list=user_model.marker_list;
MARKER_MONITORING=user_model.MARKER_MONITORING;

element_empty=user_model.element_empty;
edge_empty=HATSEdge();

% process moniter marker
% search all moniter marker nearby element
% generate edge list of all monitoring marker as hash table
% due to hash table(edge_list) property
% normal vector consist process will be done at the same time
edge_list=repmat(edge_empty,size(point_list,1),1);
for moniter_index=1:length(MARKER_MONITORING)
    [marker_elemenet,marker_index]=getMarkerElement(MARKER_MONITORING{moniter_index},marker_list);
    for element_index=1:marker_list(marker_index).element_number
        element=marker_elemenet(element_index);
        point_index_list=element.point_index_list;
        
        for point_index=1:length(point_index_list)
            % start from one edge
            vertex_index=point_index_list(point_index);
            if point_index == length(point_index_list)
                vertex_ref_index=point_index_list(1);
                % edge next
                % edge prev
            else
                vertex_ref_index=point_index_list(point_index+1);
                % edge next
                % edge prev
            end

            % add point reference of edge
            if (edge_list(vertex_index) == edge_empty)
                % means this vertex is empty
                edge=HATSEdge();
                edge.edge_number=int8(1);

                edge.vertex_ref_list=zeros(4,1,'int32');
                edge.vertex_ref_list(1)=vertex_ref_index;
                edge.element_ref_list=repmat(element_empty,4,1);
                edge.element_ref_list(1)=element;

                % edge next
                % edge prev
                % edge oppo(1)

                edge.insertRef(vertex_ref_index,1);

                edge_list(vertex_index)=edge;

                element.edge_ref_list(point_index)=1;
            else
                % means this is no empty
                edge=edge_list(vertex_index);

                % cheak ref point if exist in vertex_ref_list
                % if exist, means node index order of this element need to be rotated
                if judgeMatExistNum(edge.vertex_ref_list,vertex_ref_index)
                    % exist, rotate element node index and restart this process
                    element.point_index_list=fliplr(point_index_list);
                    element_index=element_index-1;
                    continue; % restart this process
                end

%                 % if no exist, add into vertex
%                 if length(edge.vertex_ref_list) == edge.edge_number
%                     % store array no enough, increase array
%                     edge.edge_number=edge.edge_number+1;
%                     edge.vertex_ref_list=[edge.vertex_ref_list;zeros(4,1,'int32')];
%                     edge.vertex_ref_list(edge.edge_number)=point_ref_index;
%                     edge.element_ref_list=[edge.element_ref_list;repmat(element_empty,4,1)];
%                     edge.element_ref_list(edge.edge_number)=element;
%                 else
%                     edge.edge_number=edge.edge_number+1;
%                     edge.vertex_ref_list(edge.edge_number)=point_ref_index;
%                     edge.element_ref_list(edge.edge_number)=element;
%                 end

                % if no exist, add into vertex
                edge_number=edge.edge_number+1;
                edge.edge_number=edge_number;

                edge.vertex_ref_list(edge_number)=vertex_ref_index;
                edge.element_ref_list(edge_number)=element;

                % edge next(edge_number)
                % edge prev(edge_number)
                % edge oppo(edge_number)

                edge.insertRef(vertex_ref_index,edge_number);

                element.edge_ref_list(point_index)=edge_number;
            end
        end
    end
end

user_model.edge_list=edge_list;
user_model.edge_empty=edge_empty;

% add normal_vector_list and area_list
for moniter_index=1:length(MARKER_MONITORING)
    [marker_elemenet,marker_index]=getMarkerElement(MARKER_MONITORING{moniter_index},marker_list);
    for element_index=1:marker_list(marker_index).element_number
        element=marker_elemenet(element_index);
        point_index_list=element.point_index_list;
        element.center_point=sum(point_list(point_index_list,1:user_model.dimension),1)/length(point_index_list);

        % calculate geomertry properity
        switch element.element_type
            case 3 % line
                dr12=point_list(point_index_list(2),1:2)-point_list(point_index_list(1),1:2);
                element.area=norm(dr12,2);
                element.normal_vector=[dr12(2),-dr12(1)]/element.area;
            case 5 % triangle
                dr12=point_list(point_index_list(2),1:3)-point_list(point_index_list(1),1:3);
                dr23=point_list(point_index_list(3),1:3)-point_list(point_index_list(2),1:3);

                % calculate norm_vector of element
                cross_vector=cross(dr12,dr23);
                length_cross_vector=norm(cross_vector,2);
                element.area=length_cross_vector/2;
                element.normal_vector=cross_vector/length_cross_vector;
            case 9 % quadrilateral
                dr13=point_list(point_index_list(3),1:3)-point_list(point_index_list(1),1:3);
                dr24=point_list(point_index_list(4),1:3)-point_list(point_index_list(2),1:3);

                % calculate norm_vector of element
                cross_vector=cross(dr13,dr24);
                length_cross_vector=norm(cross_vector,2);
                element.area=length_cross_vector/2;
                element.normal_vector=cross_vector/length_cross_vector;
        end
    end
end

% due to normal vector is continue, only need to judge one element normal
% if one element normal vector is incorrect, reverse all element
% volume_center_point_stagnation_point=ADtree_marker_element.center_point_list(1,:)-volume_center_point;
% normal_vector_project=volume_center_point_stagnation_point*HATS_element_list(1).normal_vector';
% if normal_vector_project < 0
%     % reverse all element
%     for element_index=1:element_number
%         HATS_element_list(element_index).normal_vector=-HATS_element_list(element_index).normal_vector;
%         HATS_element_list(element_index).point_index_list=fliplr(HATS_element_list(element_index).point_index_list);
%     end
% end
end

function getVextexList()
vertex_empty=HATSVertex();

% process moniter marker
% search all moniter marker nearby element
% generate vertex list of all marker as hash table
% due to hash table property
% normal vector consist process will be done at the same time
vertex_list=repmat(vertex_empty,size(point_list,1),1);
for moniter_index=1:length(MARKER_MONITORING)
    [marker_elemenet,marker_index]=getMarkerElement(MARKER_MONITORING{moniter_index},marker_list);
    for element_index=1:marker_list(marker_index).element_number
        element=marker_elemenet(element_index);
        point_index_list=element.point_index_list;
        for node_index=1:length(point_index_list)
            % start from one edge
            point_cur_index=point_index_list(node_index);
            if node_index < length(point_index_list)
                point_next_index=point_index_list(node_index+1);
            else
                point_next_index=point_index_list(1);
            end

            % add edge(nearby point)
            if (vertex_list(point_cur_index) == vertex_empty)
                % means this vertex is empty
                vertex=HATSVertex();
                vertex.point_next_list=zeros(4,1);
                vertex.point_next_list(1)=point_next_index;
                vertex.element_ref_list=repmat(element_empty,4,1);
                vertex.element_ref_list(1)=element;
                vertex.nearby_number=int8(1);
                vertex_list(point_cur_index)=vertex;
            else
                % means this is no empty
                vertex=vertex_list(point_cur_index);

                % cheak nearby point if exist in point nearby list
                % if exist, means node index order of this element need to be rotated
                if judgeMatExistNum(vertex.point_next_list,point_next_index)
                    % exist, rotate element node index and restart this process
                    element.point_index_list=fliplr(point_index_list);
                    element_index=element_index-1;
                    continue; % restart this process
                end

                % if no exist, add into vertex
                if length(vertex.point_next_list) == vertex.nearby_number
                    % store array no enough, increase array
                    vertex.nearby_number=vertex.nearby_number+1;
                    vertex.point_next_list=[vertex.point_next_list;zeros(4,1)];
                    vertex.point_next_list(vertex.nearby_number)=point_next_index;
                    vertex.element_ref_list=[vertex.element_ref_list;repmat(element_empty,4,1)];
                    vertex.element_ref_list(vertex.nearby_number)=element;
                else
                    vertex.nearby_number=vertex.nearby_number+1;
                    vertex.point_next_list(vertex.nearby_number)=point_next_index;
                    vertex.element_ref_list(vertex.nearby_number)=element;
                end
            end
        end
    end
end

user_model.vertex_list=vertex_list;
user_model.vertex_empty=vertex_empty;

end

function geoMakeConsistent(current_index,parent_index)
% make normal vector consistent
% check if normal_vector is consistent with parent normal_vector
% if not, reverse normal_vector
% base on public node index order to judge
% normal vector of element and node order are right-handed rule
point_index_list_current=HATS_element_list(current_index).point_index_list;
point_index_list_parent=HATS_element_list(parent_index).point_index_list;

% public point 1
public_index=1;
public_flag=0;
while public_index <= length(point_index_list_current)
    public_parent_index=1;
    while public_parent_index <= length(point_index_list_parent)
        if point_index_list_current(public_index) == point_index_list_parent(public_parent_index)
            public_flag=1;
            break;
        end
        public_parent_index=public_parent_index+1;
    end
    if public_flag
        break;
    end
    public_index=public_index+1;
end

% public point 2
public_index_sub=public_index-1;
if public_index_sub == 0
    public_index_sub=length(point_index_list_current);
end
public_parent_index_sub=public_parent_index-1;
if public_parent_index_sub == 0
    public_parent_index_sub=length(point_index_list_parent);
end

public_index_plus=public_index+1;
if public_index_plus > length(point_index_list_current)
    public_index_plus=1;
end
public_parent_index_plus=public_parent_index+1;
if public_parent_index_plus > length(point_index_list_parent)
    public_parent_index_plus=1;
end

% if edge direction do not reverse, normal should reverse
if point_index_list_current(public_index_sub) == point_index_list_parent(public_parent_index_sub) || ...
        point_index_list_current(public_index_plus) == point_index_list_parent(public_parent_index_plus)
    HATS_element_list(current_index).normal_vector=-HATS_element_list(current_index).normal_vector;
    HATS_element_list(current_index).point_index_list=fliplr(HATS_element_list(current_index).point_index_list);
end
end

function addContactElement(element_index)
% add children to element
%
point_index_list__=HATS_element_list(element_index).point_index_list;
% process symmetry
on_symmetry_flag=0;
if ~isempty(symmetry)
    on_symmetry_number=0;
    for point_nearby_index=1:point_number
        point_unit=point_index_list__(point_nearby_index);
        switch symmetry
            case 'XOY'
                dimension_index=3;
            case 'YOZ'
                dimension_index=1;
            case 'ZOX'
                dimension_index=2;
            otherwise
                error('preModel: unknown symmetry');
        end
        % check point if in symmetry face
        if abs(point_list(point_unit,dimension_index)) < geometry_torlance
            on_symmetry_number=on_symmetry_number+1;
        end
    end
    if on_symmetry_number == 2
        on_symmetry_flag=1;
    end
end

% find contact element index
contact_index_list=findContactElement(element_index,on_symmetry_flag);

% base on inflow vector judge which is inflow or outflow
% if vertical, contact node list is out flow list
if abs(HATS_element_list(element_index).normal_vector*vector_flow-1) < geometry_torlance
    out_index_list=contact_index_list;
    in_index_list=[];
else
    out_index_list=[];
    in_index_list=[];
    for contact_index__=1:length(contact_index_list)
        contact_index=contact_index_list(contact_index__);

        if ((ADtree_marker_element.center_point_list(contact_index,:)-...
                ADtree_marker_element.center_point_list(element_index,:))*vector_flow) > 0
            % outflow
            out_index_list=[out_index_list,contact_index];
        else
            % inflow
            in_index_list=[in_index_list,contact_index];
        end
    end
end

HATS_element_list(element_index).in_index_list=in_index_list;
HATS_element_list(element_index).out_index_list=out_index_list;
end

function contact_index_list=findContactElement(element_index,on_symmetry_flag)
% search nearby chilren element_index(node_index)
% if on_symmetry_flag == 1, total contact element should subtract 1
%
point_index_list__=HATS_element_list(element_index).point_index_list;

contact_index_list=[];
search_start_index_old=element_index;
search_end_index_old=element_index;
search_range_unit=search_range;

% loop to search contact element
% if number of out_index_list plus parent_node_number do not equal element edge number will contiune
% search range will increase if one search do no found enough element
while (length(contact_index_list)+on_symmetry_flag) ~= length(point_index_list__) &&...
        search_range_unit < element_number
    search_start_index=search_start_index_old-search_range_unit;
    search_start_index=max(search_start_index,1);

    search_end_index=search_end_index_old+search_range_unit;
    search_end_index=min(search_end_index,element_number);

    % search nearby element
    element_search_index_list=[search_start_index:search_start_index_old-1,search_end_index_old+1:search_end_index];
    for search_index=1:length(element_search_index_list)
        element_search_index=element_search_index_list(search_index);
        element_point_index_list=HATS_element_list(element_search_index).point_index_list;

        % if have two point same, mean is nearby element
        same_point_number=0;
        for point_index__=1:length(point_index_list__)
            for element_point_index__=1:length(element_point_index_list)
                if element_point_index_list(element_point_index__) ==...
                        point_index_list__(point_index__)
                    same_point_number=same_point_number+1;
                end
            end
        end

        if same_point_number == 2
            contact_index_list=[contact_index_list;element_search_index];
        end
    end

    % increase range
    search_range_unit=search_range_unit*2;
    search_start_index_old=search_start_index;
    search_end_index_old=search_end_index;
end
end
