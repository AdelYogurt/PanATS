classdef ADTreePoint < handle
    % Alternating Digital Tree(three dimension)
    properties
        % notice root dimension is 0
        point_list; % a pointer to store all data
        left_index_list; % left leaves list 'uint32'
        right_index_list; % right leaves list 'uint32'
        depth_list; % node depth in tree 'uint16'
        min_bou;
        max_bou;
        torlance; % default is 1e-12
    end
    methods
        function [ADT,index_list]=ADTreePoint(raw_point_list,min_bou,max_bou,torlance)
            % generate ADT by repeat_point_list
            %
            % input:
            % raw_point_list(point_number x dimension matrix), ...
            % min_bou, max_bou(default search in point_list), ...
            % torlance(min distance of point to point)
            %
            % output:
            % ADT, index_list(raw_point_list index in point_list)
            %
            if nargin > 0
                if nargin < 4
                    torlance=[];
                    if nargin < 3
                        max_bou=[];
                        if nargin < 2
                            min_bou=[];
                        end
                    end
                end
                index_list=zeros(size(raw_point_list,1),1);
                
                if isempty(torlance)
                    torlance=1e-12;
                end
                if isempty(min_bou)
                    min_bou=min(raw_point_list);
                end
                if isempty(max_bou)
                    max_bou=max(raw_point_list);
                end
                
                ADT.min_bou=min_bou;
                ADT.max_bou=max_bou;
                ADT.torlance=torlance;
                
                % initialize list
                ADT.point_list=raw_point_list;
                ADT.left_index_list=zeros(size(raw_point_list,1),1,'uint32');
                ADT.right_index_list=zeros(size(raw_point_list,1),1,'uint32');
                ADT.depth_list=zeros(size(raw_point_list,1),1,'uint16');
                
                ADT.depth_list(1)=0;
                index_list(1)=1;
                
                % add remain point
                for point_index=2:size(ADT.point_list,1)
                    point=ADT.point_list(point_index,:);
                    node_index=ADT.addExitPoint(point,point_index);
                    index_list(point_index)=node_index;
                end
            end
        end
        
        function add_node_index=addPoint(ADT,ref_point)
            % add octree node(notice have repeat protect)
            % input point will generate ADTreeNode and add into node_list
            % node_index, dimension is current node
            %
            node_low_bou=ADT.min_bou;
            node_up_bou=ADT.max_bou;
            
            node_index=1;
            dimension=0;
            depth=0;
            
            add_node_index=[];
            
            % loop to add point
            while isempty(add_node_index)
                % check current point if the same as added point
                if sum(abs(ref_point-ADT.point_list(node_index,:))) < (ADT.torlance+ADT.torlance+ADT.torlance)
                    add_node_index=node_index;
                else
                    % children node dimension
                    dimension=dimension+1;
                    if dimension > 3
                        dimension=1;
                    end
                    depth=depth+1;
                    
                    % judge node should in left or right
                    if ref_point(dimension) <= ...
                            (node_low_bou(dimension)+node_up_bou(dimension))/2
                        % left
                        if isempty(ADT.point_list(node_index,:))
                            % left children is null, add node
                            add_node_index=uint32(point_index);
                            ADT.left_index_list(node_index)=add_node_index;
                            ADT.depth_list(add_node_index).depth=depth;
                        else
                            % recursion
                            node_up_bou(dimension)=...
                                (node_low_bou(dimension)+node_up_bou(dimension))/2;
                            node_index=ADT.left_index_list(node_index);
                        end
                    else
                        % right
                        if isempty(ADT.node_list(node_index).node2_index)
                            % right children is null, add node
                            ADT.node_list=[ADT.node_list;ADTreeNode(ref_point)];
                            add_node_index=uint32(size(ADT.node_list,1));
                            ADT.node_list(node_index).node2_index=add_node_index;
                            ADT.node_list(add_node_index).depth=depth;
                        else
                            % recursion
                            node_low_bou(dimension)=...
                                (node_low_bou(dimension)+node_up_bou(dimension))/2;
                            node_index=ADT.node_list(node_index).node2_index;
                        end
                    end
                end
            end
        end
        
        function node_index_list=searchPoint(ADT,...
                low_bou,up_bou,...
                node_low_bou,node_up_bou,...
                node_index,dimension)
            % recursion to search node in range
            % search range is low_bou and up_bou
            % node_index, dimension is current node parameter
            %
            
            % cheak whether search range overlap node_bou
            if up_bou(1) > node_low_bou(1) && low_bou(1) < node_up_bou(1) &&...
                    up_bou(2) > node_low_bou(2) && low_bou(2) < node_up_bou(2) &&...
                    up_bou(3) > node_low_bou(3) && low_bou(3) < node_up_bou(3)
                % cheak if point of node is in range
                point=ADT.node_list(node_index).point;
                if (point(1) >= low_bou(1)) && (point(1) <= up_bou(1)) &&...
                        (point(2) >= low_bou(2)) && (point(2) <= up_bou(2)) &&...
                        (point(3) >= low_bou(3)) && (point(3) <= up_bou(3))
                    node_index_list=node_index;
                else
                    node_index_list=[];
                end
                
                % children node dimension
                dimension=dimension+1;
                if dimension > 3
                    dimension=1;
                end
                
                % check if exist node
                if ~isempty(ADT.node_list(node_index).node1_index)
                    % check left node(change up to middle)
                    node_low_bou_left=node_low_bou;
                    node_up_bou_left=node_up_bou;
                    node_up_bou_left(dimension)=...
                        (node_low_bou(dimension)+node_up_bou(dimension))/2;
                    % recursion
                    node_index_list=[node_index_list,searchNode(ADT,...
                        low_bou,up_bou,...
                        node_low_bou_left,node_up_bou_left,...
                        ADT.node_list(node_index).node1_index,dimension)];
                end
                
                % check if exist node
                if ~isempty(ADT.node_list(node_index).node2_index)
                    % check right node(change low to middle)
                    node_low_bou_right=node_low_bou;
                    node_low_bou_right(dimension)=...
                        (node_low_bou(dimension)+node_up_bou(dimension))/2;
                    node_up_bou_right=node_up_bou;
                    % recursion
                    node_index_list=[node_index_list,searchNode(ADT,...
                        low_bou,up_bou,...
                        node_low_bou_right,node_up_bou_right,...
                        ADT.node_list(node_index).node2_index,dimension)];
                end
            else
                node_index_list=[];
            end
        end
    end
    methods(Access = protected)
        function add_node_index=addExitPoint(ADT,ref_point,point_index)
            % add octree node(notice have repeat protect)
            % input point will generate ADTreeNode and add into node_list
            % node_index, dimension is current node
            %
            node_low_bou=ADT.min_bou;
            node_up_bou=ADT.max_bou;
            
            node_index=1; % search start from root
            dimension=0;
            depth=zeros(1,1,'uint16');
            
            add_node_index=[]; % actually node index(due to exit same node)
            
            % loop to add point
            while isempty(add_node_index)
                % check current point if the same as added point
                if sum(abs(ref_point-ADT.point_list(node_index,:))) < (ADT.torlance+ADT.torlance+ADT.torlance)
                    add_node_index=node_index;
                else
                    % search in leave node
                    % children node dimension
                    dimension=dimension+1;
                    if dimension > 3
                        dimension=1;
                    end
                    depth=depth+1;
                    
                    % judge node should in left or right
                    if ref_point(dimension) <= ...
                            (node_low_bou(dimension)+node_up_bou(dimension))/2
                        % left
                        left_node_index=ADT.left_index_list(node_index);
                        if left_node_index == 0
                            % left children is null, add node
                            add_node_index=uint32(point_index);
                            ADT.left_index_list(node_index)=add_node_index;
                            ADT.depth_list(point_index)=depth;
                        else
                            % recursion
                            node_up_bou(dimension)=...
                                (node_low_bou(dimension)+node_up_bou(dimension))/2;
                            node_index=left_node_index;
                        end
                    else
                        % right
                        right_node_index=ADT.right_index_list(node_index);
                        if right_node_index == 0
                            % right children is null, add node
                            add_node_index=uint32(point_index);
                            ADT.right_index_list(node_index)=add_node_index;
                            ADT.depth_list(point_index)=depth;
                        else
                            % recursion
                            node_low_bou(dimension)=...
                                (node_low_bou(dimension)+node_up_bou(dimension))/2;
                            node_index=right_node_index;
                        end
                    end
                end
            end
        end
    end
end
