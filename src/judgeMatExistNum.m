function [exist_flag,index]=judgeMatExistNum(mat,num)
% cheak num if exist in mat
% if exist, return 1 and place in mat, else 0
%
exist_flag=0;
for index=1:length(mat)
    if mat(index) == num
        exist_flag=1;
        return;
    end
end
index=[];
end