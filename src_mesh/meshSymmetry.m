function mesh_data=meshSymmetry(mesh_data,sym_direct,geom_torl)
% symmetry mesh data
%
if nargin < 3, geom_torl=[];end
if isempty(geom_torl), geom_torl=1e-6;end

switch sym_direct
    case 'X'
        sym_direct=1;
    case 'Y'
        sym_direct=2;
    case 'Z'
        sym_direct=3;
    otherwise
        error('meshSymmetry: unknown symmetry direction');
end

% process point_list
if isfield(mesh_data,'geometry')
    point_list=mesh_data.geometry.point_list;
    dimension=mesh_data.geometry.dimension;
    point_num=size(point_list,1);
    if dimension == 2 && sym_direct == 3
        point_list=[point_list,zeros(point_num,1)];
    end

    % search index of point on symmetry or not
    sym_index=find(abs(point_list(:,sym_direct)) < geom_torl);
    unsys_index=1:point_num;unsys_index(sym_index)=[];
    sym_num=length(sym_index);new_num=point_num-sym_num;
    
    % generate mapping list of symmetry element
    map_list=1:point_num;
    map_list(unsys_index)=(1:new_num)+point_num;

    % extend point_list to full scale
    new_point_list=point_list(unsys_index,:);
    new_point_list(:,sym_direct)=-new_point_list(:,sym_direct);
    point_list=[point_list;new_point_list];

    mesh_data.geometry.point_list=point_list;
end

% process marker
marker_name_list=fieldnames(mesh_data);
for marker_index=1:length(marker_name_list)
    marker_name=marker_name_list{marker_index};

    if strcmp(marker_name,'geometry')
        continue;
    end

    marker=mesh_data.(marker_name);

    if strcmpi(marker.type,'wgs')

    elseif strcmpi(marker.type,'stl')

    elseif strcmpi(marker.type,'MIXED3')
        elem_list=marker.element_list;
        num_list=marker.number_list;
        elem_num=length(num_list);

        % change each element idx
        node_idx=1:length(elem_list);
        node_idx(cumsum(num_list))=[];
        elem_list(node_idx)=map_list(elem_list(node_idx));

        % change element order
        data_idx=1;
        for elem_idx=1:elem_num  
            switch elem_list(data_idx)
                case 3
                    elem_list(data_idx+[1,2])=elem_list(data_idx+[2,1]);
                case 5
                    elem_list(data_idx+[1,2,3])=elem_list(data_idx+[3,2,1]);
                case 7
                    elem_list(data_idx+[1,2,3,4])=elem_list(data_idx+[4,3,2,1]);
                case 10
                    elem_list(data_idx+[1,2,3,4])=elem_list(data_idx+[4,2,3,1]);
                case 12
                    elem_list(data_idx+[1,2,3,4,5])=elem_list(data_idx+[1,4,3,2,5]);
                case 14
                    elem_list(data_idx+[1,2,3,4,5,6])=elem_list(data_idx+[3,2,1,6,5,4]);
                case 17
                    elem_list(data_idx+[1,2,3,4,5,6,7,8])=elem_list(data_idx+[1,4,3,2,5,8,7,6]);
                otherwise
            end
        end

        % add new element_list
        marker.element_list=[marker.element_list;elem_list];
        marker.number_list=[num_list,marker.number_list];
    else
        elem_list=marker.element_list;
        elem_list=map_list(elem_list);

        % change element order
        switch marker.ID
            case 3
                elem_list=elem_list(:,[2,1]);
            case 5
                elem_list=elem_list(:,[3,2,1]);
            case 7
                elem_list=elem_list(:,[4,3,2,1]);
            case 10
                elem_list=elem_list(:,[4,2,3,1]);
            case 12
                elem_list=elem_list(:,[1,4,3,2,5]);
            case 14
                elem_list=elem_list(:,[3,2,1,6,5,4]);
            case 17
                elem_list=elem_list(:,[1,4,3,2,5,8,7,6]);
            otherwise
        end

        % add new element_list
        marker.element_list=[marker.element_list;elem_list];
    end

    mesh_data.(marker_name)=marker;
end

end
