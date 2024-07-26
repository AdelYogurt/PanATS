function user_model=preModelPanel(config)
% prepare model for panel method
% calculate geometry properties of element
%
% copyright Adel 2023.03
%
config=config.data_dict;
user_model.config=config;

% initialize model
user_model.point_list=[];
user_model.marker_list=[];
user_model.element_list=[];

%% load mesh data
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
            mkr_morit=cell(length(mesh_filestr_list),1);

            for file_index=1:length(mesh_filestr_list)
                mesh_filestr=mesh_filestr_list{file_index};
                mesh_data_new=readMeshSTL(mesh_filestr,config.MESH_SCALE);
                mkr_name=fieldnames(mesh_data_new);mkr_name=mkr_name{1};
                mkr_morit{file_index}=mkr_name;
                mesh_data.(mkr_name)=mesh_data_new.(mkr_name);
            end

            % convert STL format into mesh format
            mesh_data=convertSTLToMesh(mesh_data);

            % for STL file, if MARKER_MONITORING than will analysis all file
            if ~isfield(config,'MARKER_MONITORING')
                config.MARKER_MONITORING=mkr_morit;
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
if ~isfield(config,'MARKER_MONITORING')
    error('preModelPanel: lack marker monitered');
else
    if ischar(config.MARKER_MONITORING)
        config.MARKER_MONITORING={config.MARKER_MONITORING};
    end

    % check MARKER_MONITORING if exist in mesh data
    for mkr_morit_idx=1:length(config.MARKER_MONITORING)
        mkr_morit=config.MARKER_MONITORING{mkr_morit_idx};
        exist_flag=0;
        msh_mkr_list=fieldnames(mesh_data);
        for mkr_idx=1:length(msh_mkr_list)
            if strcmp(msh_mkr_list{mkr_idx},mkr_morit) &&...
                    ~strcmp(msh_mkr_list{mkr_idx},'geometry')
                exist_flag=1;
                break;
            end
        end

        if ~exist_flag
            error('preModelPanel: marker monitered do not exist in mesh');
        end
    end
end

%% convert mesh_data to dimension, point_list, element_list and marker_list

dim=mesh_data.geometry.dimension;
pnt_list=mesh_data.geometry.point_list;

mkr_morit_list=config.MARKER_MONITORING;
mkr_num=length(mkr_morit_list);

% obtain element number of all monitoring maker
elem_num=0;
for mkr_idx=1:mkr_num
    mkr_name=mkr_morit_list{mkr_idx};
    mkr=mesh_data.(mkr_name);

    if mkr.ID == 20 % mixed element
        elem_num=elem_num+length(mkr.number_list);
    else
        elem_num=elem_num+size(mkr.element_list,1);
    end
end

% allocate memory 
if dim == 2, elem_list=nan(elem_num,2);
else, elem_list=nan(elem_num,4);end

% generate element_list and marker_list
elem_idx=1;
for mkr_idx=1:mkr_num
    mkr_name=mkr_morit_list{mkr_idx};
    mkr=mesh_data.(mkr_name);
    ID=mkr.ID;
    mkr_elem_list=mkr.element_list;
    mkr_num_list=mkr.number_list;

    % convert each mesh into Vertex
    if ID == 20 % mixed element
        mkr_elem_num=length(mkr.number_list);
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

        elem_list(mkr_elem_idx,1:mkr_num_list)=mkr_elem_list;
        elem_idx=elem_idx+mkr_elem_num;
    end

    % index list for element of marker
    mkr_list.(mkr_name)=mkr_elem_idx;
end

%% 

point_num=size(point_list,1);
elem_num=length(element_list);
% base on all element data to generate vertex list
% and base on vertex list to identify nearby element of element
% sort a pair of element index and point index into vertex list
% use hash function to decrease storage memory
% and construct vertex_list{pnt_base_idx}.pnt_ref_idx=elem_base_idx
vertex_list=repmat({zeros(1,4,'uint32');},point_num,2);
for elem_base_idx=1:elem_num
    elem_base=element_list(elem_base_idx);
    Node_idx=elem_base.Node_idx;
    node_num=elem_base.node_num;

    for node_idx=1:node_num
        % start from one edge
        pnt_base_idx=Node_idx(node_idx);
        if node_idx == node_num
            pnt_ref_idx=Node_idx(1);
            % edge next
            % edge prev
        else
            pnt_ref_idx=Node_idx(node_idx+1);
            % edge next
            % edge prev
        end

        vertex=vertex_list{pnt_base_idx,1};
        % base on hash function to insert vertex_ref_idx
        insert_idx=mod(pnt_ref_idx,4)+1;
        while insert_idx <= length(vertex) && vertex(insert_idx) ~= 0
            insert_idx=insert_idx+1;
        end

        % edge next(edge_number)
        % edge prev(edge_number)
        % edge oppo(edge_number)

        vertex_list{pnt_base_idx,1}(insert_idx)=pnt_ref_idx;
        vertex_list{pnt_base_idx,2}(insert_idx)=elem_base_idx;
    end
end

% with constructed vertex_list{pnt_base_idx}.pnt_ref_idx=elem_base_idx
% we can find out vertex_list{pnt_ref_idx}.pnt_base_idx=elem_ref_idx
% than, elem_base.Vertex_next(node_idx)=elem_ref_idx
% base on vertex list to generate Vertex_next of element
for elem_base_idx=1:elem_num
    elem_base=element_list(elem_base_idx);
    Node_idx=elem_base.Node_idx;
    node_num=elem_base.node_num;

    for node_idx=1:node_num
        % start from one edge
        pnt_base_idx=Node_idx(node_idx);
        if node_idx == node_num
            pnt_ref_idx=Node_idx(1);
            % edge next
            % edge prev
        else
            pnt_ref_idx=Node_idx(node_idx+1);
            % edge next
            % edge prev
        end

        % search reverse edge to get nearby element
        vertex=vertex_list{pnt_ref_idx,1};
        % base on hash function to insert vertex_ref_idx
        insert_idx=mod(pnt_base_idx,4)+1;
        while insert_idx < length(vertex) && vertex(insert_idx) ~= pnt_base_idx
            insert_idx=insert_idx+1;
        end
        elem_ref_idx=vertex_list{pnt_ref_idx,2}(insert_idx);

        % if equal, means that it has nearby element
        % if not, means that it may be place on symmetry place
        if vertex(insert_idx) == pnt_base_idx
            elem_base.Vertex_next(node_idx)=elem_ref_idx;
        end
    end
end

%% calculate basical geometry data of node

elem_num=length(element_list);
% calculate center_point_list, normal_vector_list and area_list
center_point_list=zeros(elem_num,3);
area_list=zeros(elem_num,1);
normal_vector_list=zeros(elem_num,3);
for elem_base_idx=1:elem_num
    elem=element_list(elem_base_idx);
    Node_idx=elem.Node_idx;
    node_num=elem.node_num;

    center_point_list(elem_base_idx,:)=sum(point_list(Node_idx,1:3),1)/node_num;

    % calculate geomertry properity
    switch elem.id
        case 3 % line
            dr12=point_list(Node_idx(2),1:2)-point_list(Node_idx(1),1:2);
            area_list(elem_base_idx,:)=norm(dr12,2);
            normal_vector_list(elem_base_idx,:)=[dr12(2),-dr12(1)]/area_list(elem_base_idx,:);
        case 5 % TRI_3
            dr12=point_list(Node_idx(2),1:3)-point_list(Node_idx(1),1:3);
            dr23=point_list(Node_idx(3),1:3)-point_list(Node_idx(2),1:3);

            % calculate norm_vector of element
            cross_vector=cross(dr12,dr23);
            length_cross_vector=norm(cross_vector,2);
            area_list(elem_base_idx,:)=length_cross_vector/2;
            normal_vector_list(elem_base_idx,:)=cross_vector/length_cross_vector;
        case 7 % QUAD_4
            dr13=point_list(Node_idx(3),1:3)-point_list(Node_idx(1),1:3);
            dr24=point_list(Node_idx(4),1:3)-point_list(Node_idx(2),1:3);

            % calculate norm_vector of element
            cross_vector=cross(dr13,dr24);
            length_cross_vector=norm(cross_vector,2);
            area_list(elem_base_idx,:)=length_cross_vector/2;
            normal_vector_list(elem_base_idx,:)=cross_vector/length_cross_vector;
    end
end

geometry.center_point_list=center_point_list;
geometry.area_list=area_list;
geometry.normal_vector_list=normal_vector_list;

user_model.geometry=geometry;
user_model.vertex_list=vertex_list;
user_model.element_list=element_list;
user_model.marker_list=mkr_list;
end
