clc;
% clear;
close all hidden;

point_nearby_list=[243,156,245,155]
point_back_list=[155,243,156,167]

index_cur=1;
point_index_list=point_nearby_list(1);
[exit_flag,index]=judgeMatExistNum...
    (point_nearby_list,point_back_list(index_cur));
while exit_flag % if exist, add into point list
    index_cur=index;
    point_index_list=[point_index_list,point_nearby_list(index_cur)];
    [exit_flag,index]=judgeMatExistNum...
        (point_nearby_list,point_back_list(index_cur));
end
% if not exist, add into point list
point_index_list=[point_index_list,point_back_list(index_cur)];
% backward search
[~,index_cur]=judgeMatExistNum...
    (point_back_list,point_nearby_list(1));
[exit_flag,index]=judgeMatExistNum...
    (point_back_list,point_nearby_list(index_cur));
while exit_flag
    index_cur=index;
    point_index_list=[point_back_list(index_cur),point_index_list];
    [exit_flag,index]=judgeMatExistNum...
        (point_back_list,point_nearby_list(index_cur));
end
point_index_list=[point_nearby_list(index_cur),point_index_list];
node_number=length(point_index_list);