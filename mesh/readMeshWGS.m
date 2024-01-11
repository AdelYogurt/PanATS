function mesh_data=readMeshWGS(mesh_filestr,scale)
% read mesh data from wgs file
%
% input:
% mesh_filestr(support .wgs file), scale(geometry zoom scale)
%
% output:
% point_list, part_list
%
% point_list is coordinate of all node
% mesh_data(single zone): mesh_data.(marker)
% marker: marker.type, marker.X, marker.Y, marker.Z
%
if nargin < 2
    scale=[];
end

% check file if exist
if exist(mesh_filestr,'file')~=2
    error('readMeshWGS: mesh file do not exist')
else
    mesh_file=fopen(mesh_filestr,'r');
end

if isempty(scale)
    scale=1;
end

% initial information
fgetl(mesh_file);

% loop each line
while ~feof(mesh_file)
    % part name
    string_list=strsplit(fgetl(mesh_file));
    if isempty(string_list{1})
        string_list=string_list(2:end);
    end
    marker_name=regexprep(string_list{1},{'[',']','"','''','-'},'');

    % part information
    string_list=strsplit(fgetl(mesh_file));
    if isempty(string_list{1})
        string_list=string_list(2:end);
    end
    information_list=str2double(string_list);
    part_ID=information_list(1);
    line_number=information_list(2);
    point_number=information_list(3);
    element_number=(line_number-1)*(point_number-1);

    part_element_number=0;

    lsymmetry=information_list(4);

    rotation=information_list(5:7);
    translation=information_list(8:10);
    scale_part=information_list(11:13);

    gsymmetry=information_list(14);

    % read data
    X=zeros(point_number,line_number);
    Y=zeros(point_number,line_number);
    Z=zeros(point_number,line_number);
    for line_index=1:line_number
        data_number=0;
        point_list=[];

        while data_number < point_number*3
            string_list=strsplit(fgetl(mesh_file));
            if isempty(string_list{1})
                string_list=string_list(2:end);
            end
            point_list=[point_list,str2double(string_list)];
            data_number=data_number+length(string_list);
        end

        point_list=reshape(point_list,3,point_number)';

        if scale ~= 1
            point_list=point_list*scale;
        end

        X(:,line_index)=point_list(:,1);
        Y(:,line_index)=point_list(:,2);
        Z(:,line_index)=point_list(:,3);
    end

    % process local symmetry
    switch lsymmetry
        case 1 % X-Z
            X=[X;flipud(X)];
            Y=[Y;-flipud(Y)];
            Z=[Z;flipud(Z)];
        case 2 % X-Y
            X=[X;flipud(X)];
            Y=[Y;flipud(Y)];
            Z=[Z;-flipud(Z)];
        case 3 % Y-Z
            X=[X;-flipud(X)];
            Y=[Y;flipud(Y)];
            Z=[Z;flipud(Z)];
    end

    % process rotation
    if rotation(2) ~= 0
        cRY=cos(rotation(2));sRY=sin(rotation(2));
        rotation_matrix=[
            cRY 0 sRY;
            0 1 0
            -sRY 0 cRY];
        X_old=X;Y_old=Y;Z_old=Z;
        X=rotation_matrix(1,1)*X_old+rotation_matrix(1,2)*Y_old+rotation_matrix(1,3)*Z_old;
        Y=rotation_matrix(2,1)*X_old+rotation_matrix(2,2)*Y_old+rotation_matrix(2,3)*Z_old;
        Z=rotation_matrix(3,1)*X_old+rotation_matrix(3,2)*Y_old+rotation_matrix(3,3)*Z_old;
    end

    if rotation(3) ~= 0
        cRZ=cos(rotation(3));sRZ=sin(rotation(3));
        rotation_matrix=[
            cRZ -sRZ 0
            sRZ cRZ 0
            0 0 1];
        X_old=X;Y_old=Y;Z_old=Z;
        X=rotation_matrix(1,1)*X_old+rotation_matrix(1,2)*Y_old+rotation_matrix(1,3)*Z_old;
        Y=rotation_matrix(2,1)*X_old+rotation_matrix(2,2)*Y_old+rotation_matrix(2,3)*Z_old;
        Z=rotation_matrix(3,1)*X_old+rotation_matrix(3,2)*Y_old+rotation_matrix(3,3)*Z_old;
    end

    if rotation(1) ~= 0
        cRX=cos(rotation(1));sRX=sin(rotation(1));
        rotation_matrix=[
            1 0 0;
            0 cRX -sRX
            0 sRX cRX];
        X_old=X;Y_old=Y;Z_old=Z;
        X=rotation_matrix(1,1)*X_old+rotation_matrix(1,2)*Y_old+rotation_matrix(1,3)*Z_old;
        Y=rotation_matrix(2,1)*X_old+rotation_matrix(2,2)*Y_old+rotation_matrix(2,3)*Z_old;
        Z=rotation_matrix(3,1)*X_old+rotation_matrix(3,2)*Y_old+rotation_matrix(3,3)*Z_old;
    end

    % process translation
    if rotation(1) ~= 0
        X=X+translation(1);
    end
    if rotation(2) ~= 0
        Y=Y+translation(2);
    end
    if rotation(3) ~= 0
        Z=Z+translation(3);
    end

    % process scale_part
    if scale_part(1) ~= 0
        X=X*scale_part(1);
    end
    if scale_part(2) ~= 0
        Y=Y*scale_part(2);
    end
    if scale_part(3) ~= 0
        Z=Z*scale_part(3);
    end

    mesh_data.(marker_name).X=X;
    mesh_data.(marker_name).Y=Y;
    mesh_data.(marker_name).Z=Z;
    mesh_data.(marker_name).type='wgs';

    % process global symmetry
    % if exist global, add mirror part 
    if gsymmetry
        switch gsymmetry
            case 1 % X-Z
                X=flipud(X);
                Y=-flipud(Y);
                Z=flipud(Z);
            case 2 % X-Y
                X=flipud(X);
                Y=flipud(Y);
                Z=-flipud(Z);
            case 3 % Y-Z
                X=-flipud(X);
                Y=flipud(Y);
                Z=flipud(Z);
        end
        
        mesh_data.([marker_name,'_SYM']).X=X;
        mesh_data.([marker_name,'_SYM']).Y=Y;
        mesh_data.([marker_name,'_SYM']).Z=Z;
        mesh_data.([marker_name,'_SYM']).type='wgs';
    end
end

fclose(mesh_file);
clear('mesh_file');

end