function writeMarker(file_name,permission,marker)
% write marker to file
% only support su2 format
% marker is cell include marker name, element number, g_Element
% write file is dat format
%
if isempty(permission)
    permission='w+';
end
if ~strcmp(file_name(end-3:end),'.dat')
    file_name=[file_name,'.dat'];
end

file_marker=fopen(file_name,permission);

% fprintf marker name
fprintf(file_marker,'MARKER_TAG= %s\n',marker{1});

% fprintf element number
fprintf(file_marker,'MARKER_ELEMS= %d\n',marker{2});

% fprintf g_Element
g_Element=marker{3};
for element_index=1:size(g_Element,1)
    element=g_Element(element_index,:);
    element_type=element(1);
    
    % element type
    fprintf(file_marker,'%d',element_type);
    
    % element point index
    switch(element_type)
        case 3
            point_number=2;
        case 5
            point_number=3;
        case 9
            point_number=4;
        case 10
            point_number=4;
        case 12
            point_number=8;
        case 13
            point_number=6;
        case 14
            point_number=5;
    end
    for point_index=1:point_number
        fprintf(file_marker,' %d',element(1+point_index));
    end
    fprintf(file_marker,'\n');
end
fclose(file_marker);
clear('file_airfoil');
end
