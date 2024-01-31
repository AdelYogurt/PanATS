clc;
% clear;
close all hidden;

% addpath([pwd,'\assist']);
% addpath([pwd,'\cfg']);
% addpath([pwd,'\data']);
% addpath([pwd,'\src']);
% addpath([pwd,'\input']);
% addpath([pwd,'\mesh']);

global user_model

%% slender

% user_model=preModelCFG('slender.cfg');
% user_model.INFORMATION=false(1);
% preModelPanel();
% 
% AOA_list=[0,3,5,7,10];
% Cl_list=zeros(1,length(AOA_list));
% Cd_list=zeros(1,length(AOA_list));
% LDratio_list=zeros(1,length(AOA_list));
% Cmy_list=zeros(1,length(AOA_list));
% for AOA_index=1:length(AOA_list)
%     AOA=AOA_list(AOA_index);
%     user_model.AOA=AOA;
%     
%     [Cl,Cd,LDratio,Cx,Cy,Cz,Cmx,Cmy,Cmz]=solveModelHypersonicInviscid();
%     [max_streamline_len]=solveModelStreamline();
%     solveModelBoundaryLayer();
%     [Cl,Cd,LDratio,Cx,Cy,Cz,Cmx,Cmy,Cmz]=solveModelHypersonicViscid();
% 
%     Cl_list(AOA_index)=Cl;
%     Cd_list(AOA_index)=Cd;
%     LDratio_list(AOA_index)=LDratio;
%     Cmy_list(AOA_index)=Cmy;
% end

% load('slender_exp.mat');
% 
% fig_hdl=figure();
% sgtitle('Slender Aerodynamic Verification');
% 
% axe_hdl_1=subplot(2,2,1);
% line(AOA_slender_exp,Cl_slender_exp,'Marker','o','Color','r','LineStyle','none');
% line(AOA_list,Cl_list,'Marker','*','Color','b');
% set(axe_hdl_1,'YLim',[-0.1,0.7],'FontName','times new roman','FontSize',12);
% axe_hdl_1.XTick=[0,3,5,7,10];
% xlabel('\alpha °');
% ylabel('{\itC}_L');
% box on;grid on;
% 
% axe_hdl_2=subplot(2,2,2);
% line(AOA_slender_exp,Cd_slender_exp,'Marker','o','Color','r','LineStyle','none');
% line(AOA_list,Cd_list,'Marker','*','Color','b');
% set(axe_hdl_2,'YLim',[0.1,0.5],'FontName','times new roman','FontSize',12);
% axe_hdl_2.XTick=[0,3,5,7,10];
% axe_hdl_2.YTick=[0.1,0.2,0.3,0.4,0.5];
% xlabel('\alpha °');
% ylabel('{\itC}_D');
% box on;grid on;
% 
% axe_hdl_3=subplot(2,2,3);
% line(AOA_slender_exp,LDratio_slender_exp,'Marker','o','Color','r','LineStyle','none');
% line(AOA_list,LDratio_list,'Marker','*','Color','b');
% set(axe_hdl_3,'YLim',[-0.2,1.8],'FontName','times new roman','FontSize',12);
% axe_hdl_3.XTick=[0,3,5,7,10];
% xlabel('\alpha °');
% ylabel('{\itC}_{L/D}');
% box on;grid on;
% 
% axe_hdl_4=subplot(2,2,4);
% line(AOA_slender_exp,Cmy_slender_exp,'Marker','o','Color','r','LineStyle','none');
% line(AOA_list,Cmy_list,'Marker','*','Color','b')
% set(axe_hdl_4,'YLim',[-0.05,0.25],'FontName','times new roman','FontSize',12);
% axe_hdl_4.XTick=[0,3,5,7,10];
% xlabel('\alpha °');
% ylabel('{\itC}_m');
% box on;grid on;
% 
% print(fig_hdl,'slender.emf','-dmeta','-r600');
% print(fig_hdl,'slender.png','-dpng','-r1200');

%% geometry

% user_model=preModelCFG('geo_test.cfg');
% preModelPanel();
% [area,volume]=solveGeometry();
