clc;
% clear all;
close all hidden;

addpath([pwd,'\src']);
addpath([pwd,'\input']);
addpath([pwd,'\mesh']);
addpath([pwd,'\assit']);
addpath([pwd,'\data']);
addpath([pwd,'\lib']);

global user_model

postModel()
output_post=user_model.output_post;

% [draw_X,draw_data]=getSYMdata(post_output.log_P_point_list);
% figure();
% plot(draw_X,log(draw_data/FREESTREAM_PRESSURE)/log(10)+1,'LineStyle','none','Marker','.')
% plot(draw_X,draw_data,'LineStyle','none','Marker','.')

load('blunt_data.mat')
[draw_coord,draw_data]=getSYMdata(output_post.Q_list);
draw_coord=draw_coord(:,2);
[draw_coord,index_list]=sort(draw_coord);
draw_coord=draw_coord-draw_coord(1);
draw_data=draw_data(index_list);

figure();
sgtitle('Blunt Heat Verification');
line(draw_coord/568.7e-3,draw_data/Q_0,'Marker','.');
line(x,Q_blunt,'LineStyle','none','Marker','o','Color','r');
legend('Estimate','EXP')
set(gca,'Xlim',[-0.1,1]);
set(gca,'Ylim',[-0.3,2.3]);
xlabel('x/L');
ylabel('q/q_0');

error_list=zeros(1,length(x));
for x_index=1:length(x)
    data=interpLinear(x(x_index),draw_coord/568.7e-3,draw_data/Q_0);
    error_list(x_index)=abs(data-Q_blunt(x_index))/Q_blunt(x_index);
end


