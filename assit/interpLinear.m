function data=interpLinear(x,x_list,data_list)
% a simple lineat interpolation
%
index=1;
while x > x_list(index+1)
    index=index+1;
end
d_data=(data_list(index+1)-data_list(index))/(x_list(index+1)-x_list(index));
data=data_list(index)+d_data*(x-x_list(index));
end
