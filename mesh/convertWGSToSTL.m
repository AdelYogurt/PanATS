function mesh_data=convertWGSToSTL(mesh_data)
% convert WSG format part_list into STL format part_list
%
% input:
% mesh_data
%
% output:
% mesh_data
%
marker_name_list=fieldnames(mesh_data);

for marker_index=1:length(marker_name_list)
    marker_name=marker_name_list{marker_index};
    if strcmp(marker_name,'geometry')
        continue;
    end

    marker=mesh_data.(marker_name);
    if ~strcmp(marker.type,'wgs')
        error('convertWGSToSTL: element_type of mesh is not LaWGS format');
    end

    X=marker.X;
    Y=marker.Y;
    Z=marker.Z;

    if any(size(X) ~= size(Y)) || any(size(X) ~= size(Z))
        error('convertWGSToSTL: size of X,Y,Z in mesh are not equal');
    end

    [point_number,line_number]=size(X);

    % convert each element to stl
    calculate_element=(point_number-1)*(line_number-1)*2;
    mesh_element_number=0;
    element_list=zeros(calculate_element*3,3);

    % write element
    for line_index=1:line_number-1
        for point_index=1:point_number-1

            x_list=X([point_index,point_index+1],[line_index,line_index+1]);
            y_list=Y([point_index,point_index+1],[line_index,line_index+1]);
            z_list=Z([point_index,point_index+1],[line_index,line_index+1]);

            % write first small element
            point1=[x_list(1),y_list(1),z_list(1)];
            point2=[x_list(3),y_list(3),z_list(3)];
            point3=[x_list(4),y_list(4),z_list(4)];
            point4=[x_list(2),y_list(2),z_list(2)];

            d12=point2-point1;
            d23=point3-point2;

            if norm(d12) < eps || norm(d23) < eps
                % small element degeneration to line, discard it
            else
                mesh_element_number=mesh_element_number+1;
                element_list(3*mesh_element_number-2,:)=point1;
                element_list(3*mesh_element_number-1,:)=point2;
                element_list(3*mesh_element_number,:)=point3;
            end

            d34=point4-point3;
            d41=point1-point4;

            if norm(d34) < eps || norm(d41) < eps
                % small element degeneration to line, discard it
            else
                mesh_element_number=mesh_element_number+1;
                element_list(3*mesh_element_number-2,:)=point3;
                element_list(3*mesh_element_number-1,:)=point4;
                element_list(3*mesh_element_number,:)=point1;
            end
        end
    end

    element_list=element_list(1:mesh_element_number*3,:);

    marker=rmfield(marker,'X');
    marker=rmfield(marker,'Y');
    marker=rmfield(marker,'Z');

    % sort element
    marker.element_list=element_list;
    marker.type='stl';

    mesh_data.(marker_name)=marker;
end

end