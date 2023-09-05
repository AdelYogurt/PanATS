function [element_list,marker_list]=preModelMesh(mesh_data,marker_name_list)
% base on mesh mesh_data to generate half edge mesh data structure
% all marker of marker_name_list will be generate a Vertex list
%
% abbreviation:
% elem: element, num: number, idx: index
%
% copyright Adel 2023.03
%
if nargin < 2
    marker_name_list=[];
end
if isempty(marker_name_list)
    marker_name_list=fieldnames(mesh_data);
end
if ischar(marker_name_list)
    marker_name_list={marker_name_list};
end
marker_index=1;
while marker_index <= length(marker_name_list)
    marker_name=marker_name_list{marker_index};
    if strcmp(marker_name,'geometry')
        marker_name_list(marker_index)=[];
    else
        marker_index=marker_index+1;
    end
end

marker_num=length(marker_name_list);

% obtain element number of all monitoring maker
elem_num=0;
for marker_idx=1:marker_num
    marker_name=marker_name_list{marker_idx};
    marker=mesh_data.(marker_name);

    if marker.ID == 20 % mixed element
        elem_num=elem_num+length(marker.number_list);
    else
        elem_num=elem_num+size(marker.element_list,1);
    end
end

% allocate memory 
element_empty=ElementVertex([],[]);
element_list=repmat(element_empty,elem_num,1);

% generate vertex list
elem_idx=1;
for marker_idx=1:marker_num
    marker_name=marker_name_list{marker_idx};
    marker=mesh_data.(marker_name);
    ID=marker.ID;
    elem_list=marker.element_list;
    num_list=marker.number_list;

    % convert each mesh into Vertex
    if ID == 20 % mixed element
        elem_num=length(marker.number_list);

        node_idx=1;
        for idx=1:elem_num
            node_num=num_list(idx);
            % new vertex
            element_list(elem_idx)=ElementVertex(elem_list(node_idx),elem_list(node_idx+(1:node_num)));
            elem_idx=elem_idx+1;
            node_idx=node_idx+node_num+1;
        end
    else
        elem_num=size(marker.element_list,1);

        for idx=1:elem_num
            % new vertex
            element_list(elem_idx)=ElementVertex(ID,elem_list(idx,:));
            elem_idx=elem_idx+1;
        end
    end

    % index list for element of marker
    Elem_idx=uint32((elem_idx-1-elem_num)+(1:elem_num));
    marker_list.(marker_name)=Elem_idx;
end

end
