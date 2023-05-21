clc;
% clear;
close all hidden;

addpath([pwd,'\src']);
addpath([pwd,'\input']);
addpath([pwd,'\mesh']);
addpath([pwd,'\lib']);
addpath([pwd,'\data']);
addpath([pwd,'\assit']);

global user_model

%% slender

user_model=preModelCFG('slender.cfg');
user_model.INFORMATION=false(1);
preModelPanel();

AOA_list=[0,3,5,7,10];
Cl_list=zeros(1,length(AOA_list));
Cd_list=zeros(1,length(AOA_list));
LDratio_list=zeros(1,length(AOA_list));
Cmy_list=zeros(1,length(AOA_list));
for AOA_index=1:length(AOA_list)
    AOA=AOA_list(AOA_index);
    user_model.AOA=AOA;
    
    [Cl,Cd,LDratio,Cx,Cy,Cz,Cmx,Cmy,Cmz]=solveModelHypersonicInviscid();
    [max_streamline_len]=solveModelStreamline();
    solveModelBoundaryLayer();
    [Cl,Cd,LDratio,Cx,Cy,Cz,Cmx,Cmy,Cmz]=solveModelHypersonicViscid();

    Cl_list(AOA_index)=Cl;
    Cd_list(AOA_index)=Cd;
    LDratio_list(AOA_index)=LDratio;
    Cmy_list(AOA_index)=Cmy;
end

load('slender_exp.mat');

figure();
sgtitle('Slender Aerodynamic Verification');

subplot(2,2,1);
line(AOA_slender_exp,Cl_slender_exp,'Marker','o','Color','r','LineStyle','none');
line(AOA_list,Cl_list,'Marker','*','Color','b');
set(gca,'YLim',[-0.1,0.7]);
xlabel('\alpha deg');
ylabel('C_L');

subplot(2,2,2);
line(AOA_slender_exp,Cd_slender_exp,'Marker','o','Color','r','LineStyle','none');
line(AOA_list,Cd_list,'Marker','*','Color','b');
set(gca,'YLim',[0.0,0.4]);
xlabel('\alpha deg');
ylabel('C_D');

subplot(2,2,3);
line(AOA_slender_exp,LDratio_slender_exp,'Marker','o','Color','r','LineStyle','none');
line(AOA_list,LDratio_list,'Marker','*','Color','b');
set(gca,'YLim',[-0.2,1.8]);
xlabel('\alpha deg');
ylabel('L/D');

subplot(2,2,4);
line(AOA_slender_exp,Cmy_slender_exp,'Marker','o','Color','r','LineStyle','none');
line(AOA_list,Cmy_list,'Marker','*','Color','b')
set(gca,'YLim',[-0.05,0.2]);
xlabel('\alpha deg');
ylabel('C_m');

%% blunt cone

% tic;
% user_model=preModelCFG('blunt_cone.cfg');
% preModelPanel();
% [Cl,Cd,LDratio,Cx,Cy,Cz,Cmx,Cmy,Cmz]=solveModelHypersonicInviscid();
% [max_streamline_len]=solveModelStreamline();
% solveModelBoundaryLayer();
% [Cl,Cd,LDratio,Cx,Cy,Cz,Cmx,Cmy,Cmz]=solveModelHypersonicViscid();
% [max_heat_flow]=solveModelHypersonicHeat();
% toc;
% 
% % displayModel('Q')
% 
% postModel();
% [draw_coord,Q_blunt_cone_HATS]=getSYMdata(user_model.output_post.Q_list);
% [draw_coord,Cp_blunt_cone_HATS]=getSYMdata(user_model.output_post.Cp_list);
% 
% load('blunt_cone_exp.mat');
% load('blunt_cone_SU2.mat');
% 
% X_blunt_cone_HATS=draw_coord(:,2);
% [X_blunt_cone_HATS,index]=sort(X_blunt_cone_HATS);
% Q_blunt_cone_HATS=Q_blunt_cone_HATS(index);
% Cp_blunt_cone_HATS=Cp_blunt_cone_HATS(index);
% 
% line(X_blunt_cone_exp,Q_blunt_cone_exp,'Marker','o','Color','k','LineStyle','none')
% line(X_blunt_cone_HATS/0.5687,Q_blunt_cone_HATS/Q_0,'Marker','.','Color','b')
% line(X_blunt_cone_SU2/0.5687,Q_blunt_cone_SU2/Q_0,'Marker','.','Color','r')
% set(gca,'YLim',[-0.5,2]);
% xlabel('x/L');
% ylabel('Q/Q_0');

% line(X_blunt_cone_HATS/0.5687,log(Cp_blunt_cone_HATS)/log(10),'Marker','.','Color','b')
% line(X_blunt_cone_SU2/0.5687,log(Cp_blunt_cone_SU2)/log(10),'Marker','.','Color','r')


%% geometry

% user_model=preModelCFG('geo_test.cfg');
% preModelPanel();
% [area,volume]=solveGeometry();