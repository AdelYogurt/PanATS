function [draw_coord,draw_data]=getSYMdata(data_list)
% read data on symmetry
%
global user_model

geom_torl=1e-6;

geometry=user_model.geometry;

dimension=geometry.dimension;
point_list=geometry.point_list;
element_list=user_model.element_list;

switch user_model.SYMMETRY
    case 'XOY'
        record_dimension=[1,2];
        identify_dimension=3;
    case 'YOZ'
        record_dimension=[2,3];
        identify_dimension=1;
    case 'ZOX'
        record_dimension=[3,1];
        identify_dimension=2;
    otherwise
        error('solveModelStreamline: nuknown SYMMETRY type');
end

% load marker point
mkr_idx=zeros(size(point_list,1));
marker_point_num=0;
for elem_idx=1:length(element_list)
    Node_idx=element_list(elem_idx).Node_idx;
    node_num=element_list(elem_idx).node_num;
    mkr_idx(marker_point_num+(1:node_num))=Node_idx;
    marker_point_num=marker_point_num+node_num;
end
mkr_idx=unique(mkr_idx(1:marker_point_num));

point_list=point_list(mkr_idx,:);
data_list=data_list(mkr_idx,:);

draw_coord=zeros(size(point_list,1),2);
draw_data=zeros(size(point_list,1),1);
sym_index=1;
for point_index=1:size(point_list,1)
   if abs(point_list(point_index,identify_dimension)) < geom_torl
       draw_coord(sym_index,:)=point_list(point_index,record_dimension);
       draw_data(sym_index)=data_list(point_index);
       sym_index=sym_index+1;
   end
end

draw_coord=draw_coord(1:(sym_index-1),:);
draw_data=draw_data(1:(sym_index-1),:);

end
