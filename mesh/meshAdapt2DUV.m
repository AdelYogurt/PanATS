function [U,V,Fval,node_list]=meshAdapt2DUV(fcn,low_bou,up_bou,torl,max_level,fval_num)
% 2D Omni-Tree
% adapt capture 2 dimemsion function value
% ensure error of linear interplation will less than torl
%
if nargin < 6
    fval_num=1;
end

% node_list which is a matrix store all node
% a node is a array, contain:
% level, idx_1-8(index of data_list), idx_c, node_index_1-4
% place:
% 3-8-4 or 3-4 or 3-8-4
% 1-5-2    6-7    6-c-7
%          1-2    1-5-2
% node:
% 1-2 or 3 or 3-4
%        1    1-2
% if children node is empty, left_index or right_index will be zero
list_add_num=50; % list will be extend only when list is not enough
node_list=zeros(list_add_num,14,'int64');
% data_list use to sort all float data include coodinate, function value
data_list=zeros(list_add_num,fval_num+2);

% add vertex of cell into data_list first
bou_1=[low_bou(1),low_bou(2)];
data_list(1,:)=[bou_1,fcn(bou_1)];
bou_2=[up_bou(1),low_bou(2)];
data_list(2,:)=[bou_2,fcn(bou_2)];
bou_3=[low_bou(1),up_bou(2)];
data_list(3,:)=[bou_3,fcn(bou_3)];
bou_4=[up_bou(1),up_bou(2)];
data_list(4,:)=[bou_4,fcn(bou_4)];

% create base root
node_list(1,:)=[0,1,2,3,4,0,0,0,0,0,0,0,0,0];

% create node tree from root
[node_num,data_num]=createNodeTree(1,4);
node_list=node_list(1:node_num,:);
data_list=data_list(1:data_num,:);

% generate U and V
[u_list,~,u_idx]=unique(data_list(:,1));
[v_list,~,v_idx]=unique(data_list(:,2));
idx_list=[u_idx,v_idx];

[U,V]=meshgrid(u_list,v_list);

% local data to Fval
Fval=nan(length(v_list),length(u_list),fval_num);

for data_index=1:data_num
    if all([idx_list(data_index,2),idx_list(data_index,1)] ~= 0)
        Fval(idx_list(data_index,2),idx_list(data_index,1),:)=data_list(data_index,3:end);
    end
end

% fit nan data
for rank_idx=1:length(v_list)
    for colume_idx=1:length(u_list)
        if isnan(Fval(rank_idx,colume_idx,1))
            Fval(rank_idx,colume_idx,:)=fcn([U(rank_idx,colume_idx),V(rank_idx,colume_idx)]);
        end
    end
end

    function [node_num,data_num]=createNodeTree(root_idx,data_num)
        % create quad tree
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
                [coord_c,coord_5,coord_6,coord_7,coord_8,...
                    fval_c,fval_5,fval_6,fval_7,fval_8]=calCell(node(2),node(3),node(4),node(5));
                [fval_pred_c,fval_pred_5,fval_pred_6,...
                    fval_pred_7,fval_pred_8]=calCellPred(node(2),node(3),node(4),node(5));

                % check u direction
                if any(abs(fval_5-fval_pred_5) > torl) || any(abs(fval_8-fval_pred_8) > torl)
                    add_u_flag=true(1);
                else
                    add_u_flag=false(1);
                end

                % check v direction
                if any(abs(fval_6-fval_pred_6) > torl) || any(abs(fval_7-fval_pred_7) > torl)
                    add_v_flag=true(1);
                else
                    add_v_flag=false(1);
                end

                % check center place
                if ~add_u_flag && ~add_v_flag && any(abs(fval_c-fval_pred_c) > torl)
                    add_u_flag=true(1);
                    add_v_flag=true(1);
                end

                if add_u_flag && add_v_flag
                    % add 5 data into data_list
                    data_new_idx=data_num+(1:5);
                    if data_num+5 > size(data_list,1)
                        data_list=[data_list;zeros(list_add_num,fval_num+2)];
                    end
                    data_list(data_new_idx,:)=[
                        coord_5,fval_5;
                        coord_6,fval_6;
                        coord_7,fval_7;
                        coord_8,fval_8;
                        coord_c,fval_c;];
                    node(6:10)=data_new_idx;
                    data_num=data_num+5;

                    % add 4 new node to node_list
                    node_new_idx=node_num+(1:4);
                    if node_num+4 > size(node_list,1)
                        node_list=[node_list;zeros(list_add_num,14)];
                    end
                    node_list(node_new_idx,:)=[...
                        node(1)+1,node(2),node(6),node(7),node(10),0,0,0,0,0,0,0,0,0;
                        node(1)+1,node(6),node(3),node(10),node(8),0,0,0,0,0,0,0,0,0;
                        node(1)+1,node(7),node(10),node(4),node(9),0,0,0,0,0,0,0,0,0;
                        node(1)+1,node(10),node(8),node(9),node(5),0,0,0,0,0,0,0,0,0;];
                    node(11:14)=node_new_idx;
                    node_num=node_num+4;

                elseif add_u_flag
                    % add 2 data into data_list
                    data_new_idx=data_num+(1:2);
                    if data_num+2 > size(data_list,1)
                        data_list=[data_list;zeros(list_add_num,fval_num+2)];
                    end
                    data_list(data_new_idx,:)=[
                        coord_5,fval_5;
                        coord_8,fval_8;];
                    node([6,9])=data_new_idx;
                    data_num=data_num+2;

                    % add 2 new node to node_list
                    node_new_idx=node_num+(1:2);
                    if node_num+2 > size(node_list,1)
                        node_list=[node_list;zeros(list_add_num,14)];
                    end
                    node_list(node_new_idx,:)=[...
                        node(1)+1,node(2),node(6),node(4),node(9),0,0,0,0,0,0,0,0,0;
                        node(1)+1,node(6),node(3),node(9),node(5),0,0,0,0,0,0,0,0,0;];
                    node([11,12])=node_new_idx;
                    node_num=node_num+2;

                elseif add_v_flag
                    % add 2 data into data_list
                    data_new_idx=data_num+(1:2);
                    if data_num+2 > size(data_list,1)
                        data_list=[data_list;zeros(list_add_num,fval_num+2)];
                    end
                    data_list(data_new_idx,:)=[
                        coord_6,fval_6;
                        coord_7,fval_7;];
                    node([7,8])=data_new_idx;
                    data_num=data_num+2;

                    % add 2 new node to node_list
                    node_new_idx=node_num+(1:2);
                    if node_num+2 > size(node_list,1)
                        node_list=[node_list;zeros(list_add_num,14)];
                    end
                    node_list(node_num+(1:2),:)=[...
                        node(1)+1,node(2),node(3),node(7),node(8),0,0,0,0,0,0,0,0,0;
                        node(1)+1,node(7),node(8),node(4),node(5),0,0,0,0,0,0,0,0,0;];
                    node([11,13])=node_new_idx;
                    node_num=node_num+2;
                else
                    node_new_idx=[];
                end

                % add to stack
                stack=[stack,node_new_idx];
                node_list(node_idx,:)=node;
            end
        end
    end

    function [coord_c,coord_5,coord_6,coord_7,coord_8,...
            fval_c,fval_5,fval_6,fval_7,fval_8]=calCell(vidx_1,vidx_2,vidx_3,vidx_4)
        % calculate index of c, 5, 6, 7, 8 place function value
        % abbreviation:
        % vidx: vertex index
        %
        coord_c=(data_list(vidx_1,[1,2])+data_list(vidx_4,[1,2]))/2;
        fval_c=fcn(coord_c);

        coord_5=(data_list(vidx_1,[1,2])+data_list(vidx_2,[1,2]))/2;
        fval_5=fcn(coord_5);
        coord_6=(data_list(vidx_1,[1,2])+data_list(vidx_3,[1,2]))/2;
        fval_6=fcn(coord_6);
        coord_7=(data_list(vidx_2,[1,2])+data_list(vidx_4,[1,2]))/2;
        fval_7=fcn(coord_7);
        coord_8=(data_list(vidx_3,[1,2])+data_list(vidx_4,[1,2]))/2;
        fval_8=fcn(coord_8);
    end

    function [fval_pred_c,fval_pred_5,fval_pred_6,...
            fval_pred_7,fval_pred_8]=calCellPred(vidx_1,vidx_2,vidx_3,vidx_4)
        % calculate index of c, 5, 6, 7, 8 place linear predict function value
        %
        fval_pred_c=(data_list(vidx_1,3:end)+data_list(vidx_2,3:end)+data_list(vidx_3,3:end)+data_list(vidx_4,3:end))/4;
        fval_pred_5=(data_list(vidx_1,3:end)+data_list(vidx_2,3:end))/2;
        fval_pred_6=(data_list(vidx_1,3:end)+data_list(vidx_3,3:end))/2;
        fval_pred_7=(data_list(vidx_2,3:end)+data_list(vidx_4,3:end))/2;
        fval_pred_8=(data_list(vidx_3,3:end)+data_list(vidx_4,3:end))/2;
    end

end
