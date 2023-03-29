classdef HATSEdge < handle
    % class of half edge list of a vertex
    % realize half edge data structure
    % edge_next↑(list)
    % edge(vertex_ref↑ element_ref↑)(list)
    % edge_prev↑(list)
    % edges of a vertex are store in self edge(vertex)
    % edge_list the same size as vertex_list(point_list)
    % so the vertex of self edge is implicit storage by edge_list
    properties
        % all list store edge of this vertex,
        % so list number size start from 4, if no enough will increase
        edge_number=int8.empty();

        % pointer to reference vertex of a halfedge
        vertex_ref_list=int32.empty(); % int32
        % pointer to reference element of a halfedge
        element_ref_list=HATSElement.empty(); % HATSElement

        % pointer to next edge of a halfedge
        edge_next_list=int32.empty(); % int32
        % pointer to previous edge of a halfedge
        edge_prev_list=int32.empty(); % int32
        % pointer to opposite edge of a halfedge, which is the index of edge in HATSEdge
        edge_oppo_list=int32.empty(); % int32

    end

    properties
        % hash map for quick search index of edge in HATSEdge
        hash_map_key=zeros(6,1,'int32'); % int32
        hash_map_data=zeros(6,1,'int8'); % int8
    end

    methods
        function insertRef(self,vertex_ref,index)
            % insert vertex_ref to hash_map
            % hash function is mod(vertex_ref,6)+1
            % clash function is linear detect
            %
            place=mod(vertex_ref,6)+1; % hash function
            if length(self.hash_map_key) < place
                self.hash_map_key(place)=vertex_ref;
                self.hash_map_data(place)=index;
            else
                while place <= length(self.hash_map_key) &&...
                        self.hash_map_key(place)
                    place=place+1; % clash function
                end
                self.hash_map_key(place)=vertex_ref;
                self.hash_map_data(place)=index;
            end
        end
        function data=getRefIndex(self,vertex_ref)
            % base on vertex ref to get index of edge in HATSEdge
            % use hash_map to realize
            %
            place=mod(vertex_ref,6)+1; % hash function
            while place < length(self.hash_map_key) &&...
                    self.hash_map_key(place) ~= vertex_ref
                place=place+1; % clash function
            end
            if self.hash_map_key(place) ~= vertex_ref
                %                 error('HATSEdge::getRefIndex: hash function error')
                data=[];
            else
                data=self.hash_map_data(place);
            end
        end
    end
end