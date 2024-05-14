function [x_list,fval_list,node_list]=meshAdapt1D(fcn,low_bou,up_bou,torl,min_level,max_level)
% Binary-tree
% adapt capture 1 dimemsion function value
% ensure error of linear interplation will less than torl
%

% node_list which is a matrix store all node
% a node is a array, contain level, index_1, index_2, index_c, node_index_1, node_index_2
% place:
% 1-c-2
% cell:
% 1-2
% if children node is empty, left_index or right_index will be zero
list_add_num=32; % node_list will be extend only when node_list is not enough
node_list=zeros(list_add_num,6,'int64');

% data_list use to sort all float data include coodinate, function value
fval_num=numel(fcn((low_bou+up_bou)/2));
data_list=zeros(list_add_num+1,fval_num+1);

% add vertex of cell into data_list first
data_list(1,:)=[low_bou,fcn(low_bou)];
data_list(2,:)=[up_bou,fcn(up_bou)];

% create base root
node_list(1,:)=[0,1,2,0,0,0];

[node_num,data_num]=createNodeTree(1,2); % create node tree from root
node_list=node_list(1:node_num,:);
data_list=data_list(1:data_num,:);

% [x_list,fval_list]=traversalInorder(1); % from small to large get list
% % add boundary info
% x_list=[data_list(1,1);x_list;data_list(2,1)];
% fval_list=[data_list(1,2:end);fval_list;data_list(2,2:end)];

x_list=data_list(:,1);
fval_list=data_list(:,2:end);

    function [node_num,data_num]=createNodeTree(root_idx,data_num)
        % create node tree from root
        %
        stack=root_idx;
        node_num=root_idx;

        while ~isempty(stack)
            % current node information
            node_idx=stack(end);
            node=node_list(node_idx,:);
            stack=stack(1:end-1);

            if node(1) < max_level
                % judge if linear predict if accptable
                % if not, create children cell
                %
                coord_c=(data_list(node(2),1)+data_list(node(3),1))/2;
                fval_c=fcn(coord_c);
                fval_c_pred=(data_list(node(2),2:end)+data_list(node(3),2:end))/2;
                
                % add 1 new data into data_list
                data_new_idx=data_num+1;
                if data_num+1 > size(data_list,1)
                    data_list=[data_list;zeros(list_add_num,fval_num+1)];
                end
                data_list(data_new_idx,:)=[coord_c,fval_c];
                node(4)=data_new_idx;
                data_num=data_num+1;

                % add 2 new node to node_list
                node_new_idx=node_num+[1,2];
                if node_num+2 > size(node_list,1)
                    node_list=[node_list;zeros(list_add_num,6)];
                end
                node_list(node_new_idx,:)=[
                    node(1)+1,node(2),node(4),0,0,0;
                    node(1)+1,node(4),node(3),0,0,0];
                node([5,6])=node_new_idx;
                node_num=node_num+2;

                if any(abs(fval_c-fval_c_pred) > torl) || node(1) < min_level-1
                    % add to stack for refine gird
                    stack=[stack,node_new_idx];
                    node_list(node_idx,:)=node;
                end
            end
        end
    end

    function [x_list,fval_list]=traversalInorder(root_idx)
        % inorder traversal node tree to obtain data
        % inorder will make sure x_list is from small to large
        %
        stack=[];
        x_list=[];
        fval_list=[];
        node_idx=root_idx;

        while ~isempty(stack) || node_idx ~= 0
            while node_idx ~= 0
                stack=[stack,node_idx];

                % node=node.left;
                node_idx=node_list(node_idx,5);
            end

            node_idx=stack(end);
            stack=stack(1:end-1);
            data_idx=node_list(node_idx,4);
            if data_idx ~= 0
                x_list=[x_list;data_list(data_idx,1)];
                fval_list=[fval_list;data_list(data_idx,2:end)];
            end

            % node=node.right;
            node_idx=node_list(node_idx,6);
        end
    end
end
