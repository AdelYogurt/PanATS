classdef ElementVertex < handle
    % class of Cell, which is the center of element
    % using half edge mesh data structure
    % notice point(self) to nearby point and element(nearby) comply with the right hand rule
    % realize half edge data structure
    % edge_next↑(list)
    % edge(vertex_ref↑ element_ref↑)(list)
    % edge_prev↑(list)
    % edges of a vertex are store in self edge(vertex)
    % edge_list the same size as vertex_list(point_list)
    % so the vertex of self edge is implicit storage by edge_list
    %
    %
    properties
        % element id of type
        id=[]; % uint8
        % rank vector, node list of element
        Node_idx=[]; % uint32
        % number of node
        node_num=[]; % double
    end
    properties
        % index to next edge of a halfedge
        Vertex_next=uint32.empty(); % uint32
        % index to previous edge of a halfedge
        Vertex_prev=uint32.empty(); % uint32
        % index to opposite edge of a halfedge
        Vertex_oppo=uint32.empty(); % uint32
    end
    methods
        function self=ElementVertex(id,Node_idx)
            self.id=id;
            self.Node_idx=Node_idx;
            node_num=(length(Node_idx));
            self.node_num=node_num;

            self.Vertex_next=zeros(1,node_num,'uint32');
            self.Vertex_prev=zeros(1,node_num,'uint32');
            self.Vertex_oppo=zeros(1,node_num,'uint32');
        end
    end
end
