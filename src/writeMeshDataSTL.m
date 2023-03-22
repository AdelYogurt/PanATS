function writeMeshDataSTL...
    (filename_mesh,data_X,data_Y,data_Z)
% write Binary STL mesh file;
% data_X, data_Y, data_Z are 1 x face_number cell
%
if length(data_X) ~= length(data_Y) || length(data_Y) ~= length(data_Z)
   error('writeMeshDataSTL: data size not equal');
end

if strcmp(filename_mesh(end-3:end),'.stl')
    filename_mesh=filename_mesh(1:end-4);
end

if ~isstruct(data_X),data_X={data_X};end
if ~isstruct(data_Y),data_Y={data_Y};end
if ~isstruct(data_Z),data_Z={data_Z};end

face_number=length(data_X);
file_mesh=fopen([filename_mesh,'.stl'],'w');

% write name to file
name_length=length(filename_mesh);
fwrite(file_mesh,['solid ',filename_mesh],'char');
fwrite(file_mesh,32*ones(80-6-name_length,1,'int8'),'int8');

% element number
fwrite(file_mesh,0,'uint32');

% write each face to file
total_element=0;
for face_index=1:face_number
    X_matrix=data_X{face_index};
    Y_matrix=data_Y{face_index};
    Z_matrix=data_Z{face_index};
    
    [rank_number,colume_number]=size(X_matrix);
    
    % check matrix data size
    if size(Y_matrix,1) ~= rank_number || size(Y_matrix,2) ~= colume_number ||...
            size(Z_matrix,1) ~= rank_number || size(Z_matrix,2) ~= colume_number
        fclose(file_mesh);
        clear('file_mesh');
        error('writeMeshDataSTL: matrix size not equal');
    end
    
    % write element
    for rank_index=1:rank_number-1
        for colume_index=1:colume_number-1
            x_list=X_matrix([rank_index,rank_index+1],[colume_index,colume_index+1]);
            y_list=Y_matrix([rank_index,rank_index+1],[colume_index,colume_index+1]);
            z_list=Z_matrix([rank_index,rank_index+1],[colume_index,colume_index+1]);
            
            % write first small element
            point1=[x_list(1),y_list(1),z_list(1)];
            point2=[x_list(2),y_list(2),z_list(2)];
            point3=[x_list(4),y_list(4),z_list(4)];
            point4=[x_list(3),y_list(3),z_list(3)];
            
            d12=point2-point1;
            d23=point3-point2;
            
            if norm(d12) < eps || norm(d23) < eps
                % small element degeneration to line, discard it
            else
                normal_vector=cross(d12,d23);
                fwrite(file_mesh,normal_vector,'float32');
                fwrite(file_mesh,point1,'float32');
                fwrite(file_mesh,point2,'float32');
                fwrite(file_mesh,point3,'float32');
                fwrite(file_mesh,[0,0],'int8');
                total_element=total_element+1;
            end
            
            d34=point4-point3;
            d41=point1-point4;
            
            if norm(d34) < eps || norm(d41) < eps
                % small element degeneration to line, discard it
            else
                normal_vector=cross(d34,d41);
                fwrite(file_mesh,normal_vector,'float32');
                fwrite(file_mesh,point3,'float32');
                fwrite(file_mesh,point4,'float32');
                fwrite(file_mesh,point1,'float32');
                fwrite(file_mesh,[0,0],'int8');
                total_element=total_element+1;
            end
        end
    end
end

% write element number
fseek(file_mesh,80,-1);
fwrite(file_mesh,total_element,'uint32');

fclose(file_mesh);
clear('file_mesh');
end