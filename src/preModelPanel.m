function preModelPanel()
% prepare model for panel method
% calculate geometry properties of element
%
% abbreviation:
% elem: element, num: number, idx: index
%
% copyright Adel 2023.03
%
global user_model

geometry=user_model.geometry;

dimension=geometry.dimension;
point_list=geometry.point_list;
element_list=user_model.element_list;

point_num=size(point_list,1);
elem_num=length(element_list);
% base on all element to generate vertex list
% and base on vertex list to identify nearby element of element
% reference vertex of vertex are sort in vertex list as hash table
vertex_list=repmat({zeros(1,4,'uint32');},point_num,2);
for elem_idx=1:elem_num
    elem=element_list(elem_idx);
    Node_idx=elem.Node_idx;

    for node_idx=1:length(Node_idx)
        % start from one edge
        vertex_base_idx=Node_idx(node_idx);
        if node_idx == length(Node_idx)
            vertex_ref_idx=Node_idx(1);
            % edge next
            % edge prev
        else
            vertex_ref_idx=Node_idx(node_idx+1);
            % edge next
            % edge prev
        end

        vertex=vertex_list{vertex_base_idx,1};
        % base on hash function to insert vertex_ref_idx
        insert_idx=mod(vertex_ref_idx,4)+1;
        while insert_idx <= length(vertex) && vertex(insert_idx) ~= 0
            insert_idx=insert_idx+1;
        end

        % edge next(edge_number)
        % edge prev(edge_number)
        % edge oppo(edge_number)

        vertex_list{vertex_base_idx,1}(insert_idx)=vertex_ref_idx;
        vertex_list{vertex_base_idx,2}(insert_idx)=elem_idx;
    end
end

% base on vertex list to generate Vertex_next of element
for elem_idx=1:elem_num
    elem=element_list(elem_idx);
    Node_idx=elem.Node_idx;

    for node_idx=1:length(Node_idx)
        % start from one edge
        vertex_base_idx=Node_idx(node_idx);
        if node_idx == length(Node_idx)
            vertex_ref_idx=Node_idx(1);
            % edge next
            % edge prev
        else
            vertex_ref_idx=Node_idx(node_idx+1);
            % edge next
            % edge prev
        end

        % search reverse edge to get nearby element
        vertex=vertex_list{vertex_ref_idx,1};
        % base on hash function to insert vertex_ref_idx
        insert_idx=mod(vertex_base_idx,4)+1;
        while insert_idx < length(vertex) && vertex(insert_idx) ~= vertex_base_idx
            insert_idx=insert_idx+1;
        end

        % if not equal, mean there have no nearby element
        if vertex(insert_idx) == vertex_base_idx
            elem.Vertex_next(node_idx)=vertex_list{vertex_ref_idx,2}(insert_idx);
        end
        
    end
end

elem_num=length(element_list);
% calculate center_point_list, normal_vector_list and area_list
center_point_list=zeros(elem_num,3);
area_list=zeros(elem_num,1);
normal_vector_list=zeros(elem_num,3);
for elem_idx=1:elem_num
    elem=element_list(elem_idx);
    Node_idx=elem.Node_idx;
    node_num=elem.node_num;

    center_point_list(elem_idx,:)=sum(point_list(Node_idx,1:3),1)/node_num;

    % calculate geomertry properity
    switch elem.id
        case 3 % line
            dr12=point_list(Node_idx(2),1:2)-point_list(Node_idx(1),1:2);
            area_list(elem_idx,:)=norm(dr12,2);
            normal_vector_list(elem_idx,:)=[dr12(2),-dr12(1)]/area_list(elem_idx,:);
        case 5 % TRI_3
            dr12=point_list(Node_idx(2),1:3)-point_list(Node_idx(1),1:3);
            dr23=point_list(Node_idx(3),1:3)-point_list(Node_idx(2),1:3);

            % calculate norm_vector of element
            cross_vector=cross(dr12,dr23);
            length_cross_vector=norm(cross_vector,2);
            area_list(elem_idx,:)=length_cross_vector/2;
            normal_vector_list(elem_idx,:)=cross_vector/length_cross_vector;
        case 6 % QUAD_4
            dr13=point_list(Node_idx(3),1:3)-point_list(Node_idx(1),1:3);
            dr24=point_list(Node_idx(4),1:3)-point_list(Node_idx(2),1:3);

            % calculate norm_vector of element
            cross_vector=cross(dr13,dr24);
            length_cross_vector=norm(cross_vector,2);
            area_list(elem_idx,:)=length_cross_vector/2;
            normal_vector_list(elem_idx,:)=cross_vector/length_cross_vector;
    end
end

geometry.center_point_list=center_point_list;
geometry.area_list=area_list;
geometry.normal_vector_list=normal_vector_list;

user_model.geometry=geometry;
end
