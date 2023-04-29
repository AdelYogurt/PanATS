function part_list=convertWGSToSTL(part_list)
% convert WSG format part_list into STL format part_list
%
if ~iscell(part_list)
    part_list={part_list};
end

for part_index=1:length(part_list)
    part=part_list{part_index};
    mesh_list=part.mesh_list;
    part_element_number=0;

    for mesh_index=1:length(mesh_list)
        mesh=mesh_list{mesh_index};
        if ~strcmp(mesh.element_type,'wgs')
            error('convertWGSToSTL: element_type of mesh is not LaWGS format');
        end

        X=mesh.X;
        Y=mesh.Y;
        Z=mesh.Z;
        
        if any(size(X) ~= size(Y)) || any(size(X) ~= size(Z))
            error('convertWGSToSTL: size of X,Y,Z in mesh are not equal');
        end

        [point_number,line_number]=size(X);

        % convert each element to stl
        calculate_element=(point_number-1)*(line_number-1)*2;
        element_number=0;
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
                    element_number=element_number+1;
                    element_list(3*element_number-2,:)=point1;
                    element_list(3*element_number-1,:)=point2;
                    element_list(3*element_number,:)=point3;
                end

                d34=point4-point3;
                d41=point1-point4;

                if norm(d34) < eps || norm(d41) < eps
                    % small element degeneration to line, discard it
                else
                    element_number=element_number+1;
                    element_list(3*element_number-2,:)=point3;
                    element_list(3*element_number-1,:)=point4;
                    element_list(3*element_number,:)=point1;
                end
            end
        end

        element_list=element_list(1:element_number*3,:);

        mesh=rmfield(mesh,'X');
        mesh=rmfield(mesh,'Y');
        mesh=rmfield(mesh,'Z');

        % sort element
        mesh.element_list=element_list;
        mesh.element_type='stl';
        part_element_number=part_element_number+element_number;

        mesh_list{mesh_index}=mesh;
    end

    part.mesh_list=mesh_list;
    part_list{part_index}=part;
    part_list{part_index}.element_number=part_element_number;
end

if length(part_list) == 1
    part_list=part_list{1};
end

end