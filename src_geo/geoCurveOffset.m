function [curve,line_flag]=geoCurveOffset(curve,offset)
% shifting line with offset
% direction is outside
%
origin=mean(curve,1);
curve=curve-origin;
line_flag=false;
edge_number=size(curve,1);

% calculate line data and move line
edge_data_list=zeros(edge_number,3); % A, B, C
for edge_index=1:edge_number
    point_1=curve(edge_index,:);
    if edge_index == size(curve,1)
        point_2=curve(1,:);
    else
        point_2=curve(edge_index+1,:);
    end

    dr=point_2-point_1;
    
    edge_data_list(edge_index,1)=-dr(2); % A (dr_y=-A)
    edge_data_list(edge_index,2)=dr(1); % B (dr_x=B)

    if norm(dr) == 0
        return;
    end

    % move line
    edge_data_list(edge_index,3)=point_1(1)*point_2(2)-point_1(2)*point_2(1)+...
        offset/norm(dr)*sum(dr.^2); % C
end

% solve new cross point
for edge_index=1:edge_number
    if edge_index == 1
        edge_prev_index=edge_number;
    else
        edge_prev_index=edge_index-1;
    end

    matrix=[edge_data_list(edge_index,1),edge_data_list(edge_index,2);
        edge_data_list(edge_prev_index,1),edge_data_list(edge_prev_index,2)];
    if det(matrix) < eps
        line_flag=true;
        curve=[];
        return;
    end

    curve(edge_index,:)=matrix\[-edge_data_list(edge_index,3);-edge_data_list(edge_prev_index,3)];
end

curve=curve+origin;
end
