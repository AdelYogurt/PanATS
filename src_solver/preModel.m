function [model,geom_out]=preModel(model)
% prepare model for panel method
% load mesh data from file
% create connectivity of element
% calculate geometry properties of element
%
% copyright Adel 2023.03
%
config=model.config; % load config
geometry=model.geometry; % load geometry

INFORMATION=config.INFORMATION; % whether print information

% load config
SYMMETRY=config.SYMMETRY;

% initialize geometry
geometry.point_list=[];
geometry.marker_list=[];
geometry.element_list=[];
geometry.connectivity_matrix=[];
geometry.normal_vector_list=[];
geometry.area_list=[];

geom_tol=1e-6;
switch SYMMETRY
    case 'XOY'
        identify_dim=3;
    case 'YOZ'
        identify_dim=1;
    case 'ZOX'
        identify_dim=2;
    otherwise
        identify_dim=0;
end

%% load mesh data

if INFORMATION
    fprintf('Loading mesh form %s.\n',config.MESH_FILENAME);
    fprintf('Loading mesh format is %s.\n',config.MESH_FORMAT);
end

if ~isfield(config,'mesh_data') || isempty(config.mesh_data)
    switch config.MESH_FORMAT
        case 'SU2'
            % only need marker infomation
            mesh_data=readMeshSU2(config.MESH_FILENAME,config.MESH_SCALE,false(1),true(1));
        case 'STL'
            mesh_filestr_list=config.MESH_FILENAME;
            if ischar(mesh_filestr_list)
                mesh_filestr_list={mesh_filestr_list};
            end
            mesh_data=struct();
            mkr_name_list=cell(length(mesh_filestr_list),1);

            for file_index=1:length(mesh_filestr_list)
                mesh_filestr=mesh_filestr_list{file_index};
                mesh_data_new=readMeshSTL(mesh_filestr,config.MESH_SCALE);
                mkr_name=fieldnames(mesh_data_new);mkr_name=mkr_name{1};
                mkr_name_list{file_index}=mkr_name;
                mesh_data.(mkr_name)=mesh_data_new.(mkr_name);
            end

            % convert STL format into mesh format
            mesh_data=convertSTLToMesh(mesh_data);

            % for STL file, if MARKER_MONITORING than will analysis all file
            if isempty(config.MARKER_MONITORING)
                config.MARKER_MONITORING=mkr_name_list;
            end
        case 'INP'
            % read SU2 format mesh data into mesh format
            mesh_data=readMeshINP(config.MESH_FILENAME,config.MESH_SCALE);
        case 'WGS'
            mesh_data=readMeshWGS(config.MESH_FILENAME,config.MESH_SCALE);

            % convert LaWGS format into mesh format
            mesh_data=convertWGSToMesh(mesh_data);
        case 'CGNS'
            mesh_data=readMeshCGNS(config.MESH_FILENAME,config.MESH_SCALE);
    end
else
    mesh_data=config.mesh_data;
    switch config.MESH_FORMAT
        case 'WGS'
            % convert LaWGS format into mesh format
            mesh_data=convertWGSToMesh(mesh_data);
        case 'STL'
            % convert STL format into mesh format
            mesh_data=convertSTLToMesh(mesh_data);
    end
end

% process marker monitoring
if ~isempty(config.MARKER_MONITORING) && ischar(config.MARKER_MONITORING)
    config.MARKER_MONITORING={config.MARKER_MONITORING};
end

% check MARKER_MONITORING if exist in mesh data
for mkr_morit_idx=1:length(config.MARKER_MONITORING)
    mkr_monit=config.MARKER_MONITORING{mkr_morit_idx};

    exist_flag=0;
    msh_mkr_list=fieldnames(mesh_data);
    for mkr_idx=1:length(msh_mkr_list)
        if strcmp(msh_mkr_list{mkr_idx},mkr_monit) &&...
                ~strcmp(msh_mkr_list{mkr_idx},'geometry')
            exist_flag=1;
            break;
        end
    end

    if ~exist_flag
        error('preModel: marker monitered do not exist in mesh');
    end
end

if INFORMATION
    if ~isempty(config.MARKER_MONITORING)
        fprintf('Marker monitering is/are ( %s ).\n',strjoin(config.MARKER_MONITORING,', '));
    end
end

%% convert mesh_data to dimension, point_list, element_list and marker_list

dim=mesh_data.geometry.dimension;
pnt_list=mesh_data.geometry.point_list;
mesh_data=rmfield(mesh_data,'geometry');

len_x=max(pnt_list(:,1))-min(pnt_list(:,1));
len_y=max(pnt_list(:,2))-min(pnt_list(:,2));
if dim == 3,len_z=max(pnt_list(:,3))-min(pnt_list(:,3));end

mkr_name_list=fieldnames(mesh_data);
mkr_num=length(mkr_name_list);

% load all element for mesh_data
elem_num=0;
for mkr_idx=1:mkr_num
    mkr_name=mkr_name_list{mkr_idx};
    mkr=mesh_data.(mkr_name);

    % for different dimension load different element
    if dim == 2 && mkr.ID ~= 3,continue;end
    if dim == 3 && mkr.ID ~= 5 && mkr.ID ~= 7 && mkr.ID ~= 20,continue;end

    if mkr.ID == 20 % mixed element
        elem_num=elem_num+length(mkr.number_list);
    else
        elem_num=elem_num+size(mkr.element_list,1);
    end
end

% allocate memory 
if dim == 2, elem_list=nan(elem_num,2);
else, elem_list=nan(elem_num,4);end
mkr_list=struct();

% generate element_list and marker_list
elem_idx=1;
for mkr_idx=1:mkr_num
    mkr_name=mkr_name_list{mkr_idx};
    mkr=mesh_data.(mkr_name);

    % for different dimension load different element
    if dim == 2 && mkr.ID ~= 3,continue;end
    if dim == 3 && mkr.ID ~= 5 && mkr.ID ~= 7 && mkr.ID ~= 20,continue;end

    % convert each mesh into Vertex
    ID=mkr.ID;
    mkr_elem_list=mkr.element_list;
    if ID == 20 % mixed element
        mkr_num_list=mkr.number_list;
        mkr_elem_num=length(mkr_num_list);
        mkr_elem_idx=uint32((elem_idx):(elem_idx-1+mkr_elem_num))';

        node_idx=1;
        for idx=1:mkr_elem_num
            node_num=mkr_num_list(idx);
            elem_list(elem_idx,1:node_num)=mkr_elem_list(node_idx+(1:node_num));
            elem_idx=elem_idx+1;
            node_idx=node_idx+node_num+1;
        end
    else
        mkr_elem_num=size(mkr.element_list,1);
        mkr_elem_idx=uint32((elem_idx):(elem_idx-1+mkr_elem_num))';

        elem_list(mkr_elem_idx,1:size(mkr_elem_list,2))=mkr_elem_list;
        elem_idx=elem_idx+mkr_elem_num;
    end

    % index list for element of marker
    mkr_list.(mkr_name)=mkr_elem_idx;
end

if INFORMATION,fprintf('Model dimension is %d.\n',dim);end
clear('mesh_data');

%% create half edge data structure to sort connectivity

elem_num=size(elem_list,1);
if dim == 2
    % create all vertex
    vtx_i=reshape(elem_list,[],1);
    vtx_j=[2*ones(elem_num,1);ones(elem_num,1)];
    vtx_v=repmat(1:elem_num,2,1)';

    % create connectivity matrix
    conn_mat=sparse(vtx_i,vtx_j,vtx_v);

    % check element topology

else
    % separate triangular and quadrilateral element
    bool=isnan(elem_list(:,end));
    tri_idx=find(bool);
    quad_idx=find(~bool);

    % create all vertex on edge, direction is vtx_i to vtx_j
    vtx_i=reshape(elem_list',[],1);
    vtx_j=elem_list;
    vtx_j(tri_idx,:)=vtx_j(tri_idx,[2,3,1,4]);
    vtx_j(quad_idx,:)=vtx_j(quad_idx,[2,3,4,1]);
    vtx_j=reshape(vtx_j',[],1);
    vtx_v=repmat(1:elem_num,4,1)';
    vtx_v(tri_idx,4)=nan;
    vtx_v=reshape(vtx_v',[],1);

    nan_idx=tri_idx*4;

    vtx_i(nan_idx)=[];
    vtx_j(nan_idx)=[];
    vtx_v(nan_idx)=[];

    % create connectivity matrix
    conn_mat=sparse(vtx_i,vtx_j,vtx_v);

    % check element topology

end

%% calculate basical geometry data of element

% calculate center_point_list, point_normal_vector_list, element_normal_vector_list and area_list
pnt_num=size(pnt_list,1);
P_nmvctr_list=zeros(pnt_num,dim);
sum_weight_list=zeros(pnt_num,1);
if dim == 2 % BAR_2
    node_idx_1=elem_list(:,1);
    node_idx_2=elem_list(:,2);

    cntr_pnt_list=(pnt_list(node_idx_1,1:2)+pnt_list(node_idx_2,1:2))/2;

    % calculate normal vector of element and point
    tgvctr_list=pnt_list(node_idx_2,1:2)-pnt_list(node_idx_1,1:2);
    E_nmvctr_list=[-tgvctr_list(:,2),tgvctr_list(:,1)];
    area_list=sqrt(sum(E_nmvctr_list.^2,2));

    % add point normal vector
    P_nmvctr_list(node_idx_1,:)=P_nmvctr_list(node_idx_1,:)+E_nmvctr_list;
    sum_weight_list(node_idx_1)=sum_weight_list(node_idx_1)+area_list;
    P_nmvctr_list(node_idx_2,:)=P_nmvctr_list(node_idx_2,:)+E_nmvctr_list;
    sum_weight_list(node_idx_2)=sum_weight_list(node_idx_2)+area_list;

    E_nmvctr_list=E_nmvctr_list./area_list;
    P_nmvctr_list=P_nmvctr_list./sum_weight_list;
else % TRI_3 and QUAD_4
    node_idx_1=elem_list(:,1);
    node_idx_2=elem_list(:,2);
    node_idx_3=elem_list(:,3);
    node_idx_4=elem_list(:,4);

    cntr_pnt_list=pnt_list(node_idx_1,1:3)+pnt_list(node_idx_2,1:3)+pnt_list(node_idx_3,1:3);
    cntr_pnt_list(tri_idx,:)=cntr_pnt_list(tri_idx,:)/3;
    cntr_pnt_list(quad_idx,:)=(cntr_pnt_list(quad_idx,:)+pnt_list(node_idx_4(quad_idx),1:3))/4;

    % process vector and edge of triangular and quadrilateral element respectively

    % edges of element
    if INFORMATION,fprintf('Checking mesh quality.\n');end
    edg_len_list=inf(elem_num,4);

    edg_3_list=nan(elem_num,3);
    edg_4_list=nan(elem_num,3);
    edg_1_list=pnt_list(node_idx_2,1:3)-pnt_list(node_idx_1,1:3);
    edg_2_list=pnt_list(node_idx_3,1:3)-pnt_list(node_idx_2,1:3);
    edg_3_list(tri_idx,:)=pnt_list(node_idx_1(tri_idx,:),1:3)-pnt_list(node_idx_3(tri_idx,:),1:3);
    edg_3_list(quad_idx,:)=pnt_list(node_idx_4(quad_idx,:),1:3)-pnt_list(node_idx_3(quad_idx,:),1:3);
    edg_4_list(quad_idx,:)=pnt_list(node_idx_1(quad_idx,:),1:3)-pnt_list(node_idx_4(quad_idx,:),1:3);

    % for TRI_3
    vctr_i=edg_1_list;
    vctr_j=edg_2_list;

    % for QUAD_4
    vctr_i(quad_idx,:)=pnt_list(node_idx_3(quad_idx),1:3)-pnt_list(node_idx_1(quad_idx),1:3);
    vctr_j(quad_idx,:)=pnt_list(node_idx_4(quad_idx),1:3)-pnt_list(node_idx_2(quad_idx),1:3);

    edg_len_list(:,1)=sqrt(sum(edg_1_list.^2,2));
    edg_len_list(:,2)=sqrt(sum(edg_2_list.^2,2));
    edg_len_list(:,3)=sqrt(sum(edg_3_list.^2,2));
    edg_len_list(quad_idx,4)=sqrt(sum(edg_4_list(quad_idx,:).^2,2));

    min_len_list=min(edg_len_list,[],2);
    edg_len_list(tri_idx,4)=0;
    max_len_list=max(edg_len_list,[],2);

    AR_list=min_len_list./max_len_list;
    min_AR=min(AR_list);
    avg_AR=mean(AR_list);
    max_AR=max(AR_list);

    if INFORMATION
        fprintf('Aspect Ratio: min %f, average %f, max %f.\n',min_AR,avg_AR,max_AR);
    end

    if avg_AR > 0.5
        FLAG_AREA=true;
        if INFORMATION,fprintf('Using area contribute method to calculate normal vector of point.\n');end
    else
        FLAG_AREA=false;
        if INFORMATION,fprintf('Using angle contribute method to calculate normal vector of point.\n');end
    end

    % calculate normal vector of element and point
    E_nmvctr_list=cross(vctr_i,vctr_j,2);
    E_nmvctr_list(tri_idx,:)=E_nmvctr_list(tri_idx,:)/2; % for triangular element, area is half of cross vector length
    area_list=sqrt(sum(E_nmvctr_list.^2,2));

    if FLAG_AREA % using area contribute method
        % add point normal vector
        for elem_idx=1:elem_num
            if isnan(elem_list(elem_idx,4)),node_idx_list=elem_list(elem_idx,1:3);
            else,node_idx_list=elem_list(elem_idx,1:4);end
            P_nmvctr_list(node_idx_list,:)=P_nmvctr_list(node_idx_list,:)+E_nmvctr_list(elem_idx,:);
        end
    end

    % re-normalize normal vector of element
    E_nmvctr_list=E_nmvctr_list./area_list;
    
    if ~FLAG_AREA % using angle contribute method
        % calculate angle of element
        ang_list=zeros(elem_num,4);

        edg_1_list=edg_1_list./edg_len_list(:,1);
        edg_2_list=edg_2_list./edg_len_list(:,2);
        edg_3_list=edg_3_list./edg_len_list(:,3);
        edg_4_list(quad_idx,:)=edg_4_list(quad_idx,:)./edg_len_list(quad_idx,4);

        ang_list(tri_idx,1)=acos(sum(-edg_3_list(tri_idx,:).*edg_1_list(tri_idx,:),2));
        ang_list(quad_idx,1)=acos(sum(-edg_4_list(quad_idx,:).*edg_1_list(quad_idx,:),2));
        ang_list(:,2)=acos(sum(-edg_1_list.*edg_2_list,2));
        ang_list(:,3)=acos(sum(-edg_2_list.*edg_3_list,2));
        ang_list(quad_idx,4)=acos(sum(-edg_3_list(quad_idx,:).*edg_4_list(quad_idx,:),2));
        
        % add point normal vector
        for elem_idx=1:elem_num
            if isnan(elem_list(elem_idx,4))
                node_idx_list=elem_list(elem_idx,1:3);
                E_nmvctr_add=E_nmvctr_list(elem_idx,:).*ang_list(elem_idx,1:3)';
            else
                node_idx_list=elem_list(elem_idx,1:4);
                E_nmvctr_add=E_nmvctr_list(elem_idx,:).*ang_list(elem_idx,1:4)';
            end
            P_nmvctr_list(node_idx_list,:)=P_nmvctr_list(node_idx_list,:)+E_nmvctr_add;
        end
    end

    if identify_dim % identify which point of element on symmetry face and process symmetry
        Bool_sym=abs(pnt_list(:,identify_dim)) <= geom_tol;
        P_nmvctr_list(Bool_sym,identify_dim)=0;
    end

    % re-normalize normal vector of point
    P_nmvctr_list=P_nmvctr_list./sqrt(sum(P_nmvctr_list.^2,2));
end

%% calculate basical geometry of model

% area is sum of all element area
area=sum(area_list);
area_x=sum(area_list.*max(E_nmvctr_list(:,1),0));
area_y=sum(area_list.*max(E_nmvctr_list(:,2),0));
if dim == 3,area_z=sum(area_list.*max(E_nmvctr_list(:,3),0));end

if dim == 2 % BAR_2
    TPL=[pnt_list(elem_list(:,1),1:2),pnt_list(elem_list(:,2),1:2)];
    volume=sum((TPL(:,2).*TPL(:,3)-TPL(:,1).*TPL(:,4)))/2;
else % TRI_3 and QUAD_4
    % convert all element to triangular element, load all node list
    elem_num=length(elem_list);
    quad_num=length(quad_idx);
    TPL=zeros(elem_num+quad_num,9);
    TPL(1:elem_num,:)=[pnt_list(elem_list(:,1),1:3),pnt_list(elem_list(:,2),1:3),pnt_list(elem_list(:,3),1:3)];
    TPL((elem_num+1):(elem_num+quad_num),:)=[pnt_list(elem_list(quad_idx,1),1:3),pnt_list(elem_list(quad_idx,3),1:3),pnt_list(elem_list(quad_idx,4),1:3)];

    % volume is sum of all directed volume
    volume=sum((-TPL(:,7).*TPL(:,5).*TPL(:,3)+TPL(:,4).*TPL(:,8).*TPL(:,3)+TPL(:,7).*TPL(:,2).*TPL(:,6)+...
        -TPL(:,1).*TPL(:,8).*TPL(:,6)-TPL(:,4).*TPL(:,2).*TPL(:,9)+TPL(:,1).*TPL(:,5).*TPL(:,9)))/6;
end

%% sort data

if model.config.REF_LENGTH == 0, model.config.REF_LENGTH=len_x;end
if model.config.REF_AREA == 0, model.config.REF_AREA=area_z;end

geom_out.len_x=len_x;
geom_out.len_y=len_y;
geom_out.area=area;
geom_out.area_x=area_x;
geom_out.area_y=area_y;
geom_out.volume=volume;
if dim == 3
    geom_out.len_z=len_z;
    geom_out.area_z=area_z;
end

geometry.dimension=dim;
geometry.point_list=pnt_list;
geometry.marker_list=mkr_list;
geometry.element_list=elem_list;
geometry.connectivity_matrix=conn_mat;
geometry.center_point_list=cntr_pnt_list;
geometry.E_Nvector_list=E_nmvctr_list;
geometry.P_Nvector_list=P_nmvctr_list;
geometry.area_list=area_list;

model.geometry=geometry;
end
