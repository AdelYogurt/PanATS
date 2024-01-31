clc;
% clear;
close all hidden;

addpath([pwd,'\assist']);
addpath([pwd,'\cfg']);
addpath([pwd,'\data']);
addpath([pwd,'\input']);
addpath([pwd,'\mesh']);
addpath([pwd,'\src']);

global user_model

%% blunt cone

user_model=preModelCFG('blunt_cone.cfg');
preModelPanel();
[Cl,Cd,LDratio,Cx,Cy,Cz,Cmx,Cmy,Cmz]=solveModelHypersonicInviscid();
[max_streamline_len]=solveModelStreamline();
solveModelBoundaryLayer();
[Cl,Cd,LDratio,Cx,Cy,Cz,Cmx,Cmy,Cmz]=solveModelHypersonicViscid();
[max_heat_flow]=solveModelHypersonicHeat();
displayModel('HF')

%% post process

postModel()
output_post=user_model.output_post;

load('blunt_cone_exp.mat')
[draw_coord,draw_data]=getSYMdata(output_post.HF_list);
draw_coord=draw_coord(:,2)/568.7e-3; % scale to length
[draw_coord,index_list]=sort(draw_coord);
draw_coord=draw_coord-draw_coord(1);
draw_data=draw_data(index_list);

figure();
sgtitle('Blunt Heat Verification');
line(draw_coord,draw_data/HF_0,'Marker','.');
line(X_blunt_cone_exp,HF_blunt_cone_exp,'LineStyle','none','Marker','o','Color','r');
legend('Estimate','EXP')
set(gca,'Xlim',[-0.1,1.0]);
set(gca,'Ylim',[-0.3,2.8]);
xlabel('x/L');
ylabel('q/q_0');

interp_HF=interpLinear(draw_coord,draw_data/HF_0,X_blunt_cone_exp)';
err_list=abs(interp_HF-HF_blunt_cone_exp)./HF_blunt_cone_exp;

