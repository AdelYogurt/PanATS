function writeNode(file_name,permission,point_index_list)
% function to print point coordinate to file
% include point index(global) point coordinate
%
global g_Point
if isempty(permission)
    permission='w+';
end

if ~strcmp(file_name(end-3:end),'.dat')
    file_name=[file_name,'.dat'];
end

dimension=size(g_Point,2);

file_point=fopen(file_name,permission);

% fprintf g_Point of point_index_list
for point_index=1:size(point_index_list,1)
    point_unit=point_index_list(point_index);
    
    % point index(global)
    fprintf(file_point,'%d',point_unit);
    
    % point coordinate
    for dimension_index=1:dimension
        fprintf(file_point,' % 9.8e',g_Point(point_unit+1,dimension_index));
    end
    
    fprintf(file_point,'\n');
end
fclose(file_point);
clear('file_point');
end
