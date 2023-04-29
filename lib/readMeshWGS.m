function part_list=readMeshWGS(filename_mesh,scale,INFORMATION)
% read mash data from wgs file
% input:
% filename_mesh(support .wgs file), scale(geometry zoom scale), ...
% INFORMATION(true or false)
%
% output:
% point_list, part_list
%
% point_list is coordinate of all node
% part_list{part.name, part.mesh_list{mesh.element_type, mesh.X, mesh.Y, mesh.Z}}
%
if nargin < 3
    INFORMATION=true(1);
    if nargin < 2
        scale=[];
    end
end

% cheak filename
if length(filename_mesh) > 4
    if ~strcmpi(filename_mesh((end-3):end),'.wgs')
        filename_mesh=[filename_mesh,'.wgs'];
    end
else
    filename_mesh=[filename_mesh,'.wgs'];
end

% check file if exist
if exist(filename_mesh,'file')~=2
    error('readMeshWGS: mesh file do not exist')
else
    file_mesh=fopen(filename_mesh,'r');
end

point_list=[];
part_list={};
part_number=0;

if isempty(scale)
    scale=1;
end

if INFORMATION
    disp('readMeshWGS: read mash data begin');
end

% initial information
fgetl(file_mesh);

read_part_flag=0;

% loop each line
part_index=0;
while ~feof(file_mesh)
    % part name
    string_list=strsplit(fgetl(file_mesh));
    if isempty(string_list{1})
        string_list=string_list(2:end);
    end
    part_name=string_list{1};

    % part information
    string_list=strsplit(fgetl(file_mesh));
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
            string_list=strsplit(fgetl(file_mesh));
            if isempty(string_list{1})
                string_list=string_list(2:end);
            end
            point_list=[point_list,str2double(string_list)];
            data_number=data_number+length(string_list);
        end

        point_list=reshape(point_list,3,point_number)'*scale;

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

    mesh.X=X;
    mesh.Y=Y;
    mesh.Z=Z;
    mesh.element_type='wgs';
    mesh.element_number=element_number;
    
    mesh_list(1,1)={mesh};
    part_element_number=part_element_number+element_number;

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
        mesh.X=X;
        mesh.Y=Y;
        mesh.Z=Z;
        mesh.element_type='wgs';

        mesh_list{length(mesh_list)+1,1}=mesh;
        part_element_number=part_element_number+element_number;
    end

    part.name=part_name;
    part.mesh_list=mesh_list;
    part.element_number=part_element_number;

    part_number=part_number+1;
    part_list{part_number,1}=part;
end

fclose(file_mesh);
clear('file_mesh');

if INFORMATION
    disp('readMeshWGS: read mash data done');
end

end