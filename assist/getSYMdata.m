function [draw_coord,draw_data]=getSYMdata(data_list)
% read data on symmetry
%
global user_model

geometry_torlance=1e-6;

dimension=user_model.geometry.dimension;

point_list=user_model.point_list;
edge_list=user_model.edge_list;
marker_list=user_model.marker_list;

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

output_post=user_model.output_post;
marker_boolean=output_post.marker_boolean;
marker_point_number=sum(marker_boolean);

draw_coord=zeros(marker_point_number,2);
draw_data=zeros(marker_point_number,1);
marker_index=1;
for point_index=1:size(point_list,1)
   if marker_boolean(point_index) && abs(point_list(point_index,identify_dimension)) < geometry_torlance
       draw_coord(marker_index,:)=point_list(point_index,record_dimension);
       draw_data(marker_index)=data_list(point_index);
       marker_index=marker_index+1;
   end
end

draw_coord=draw_coord(1:(marker_index-1),:);
draw_data=draw_data(1:(marker_index-1),:);

end
