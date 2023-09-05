clc;
clear;
close all hidden;

% load('slender_exp.mat');
% load('slender_SU2.mat')
% load('slender_HATS.mat')
% 
% Cl_slender_SU2=[0,0.082,0.165,0.315,0.575];
% fig_hdl=figure(1);
% line(AOA_slender_exp,Cl_slender_exp,'Marker','s','Color','r','LineStyle','none','MarkerFaceColor','r');
% line(AOA_slender_SU2,Cl_slender_SU2,'Marker','o','Color','b','LineWidth',1,'LineStyle','-');
% line(AOA_list,Cl_list,'Marker','^','Color',[0.9290 0.6940 0.1250],'LineWidth',1,'LineStyle','--');
% set(gca,'YLim',[-0.1,0.7]);
% xlabel('\alpha/°');
% ylabel('C_L');
% legend('Exp','HF','LF','Location','northwest')
% fig_hdl.set('Position',[488,342,280,210])
% mesh_data on;
% 
% LDratio_slender_exp=[0,0.33,0.68,1.00,1.45];
% LDratio_slender_SU2=[0,0.35,0.72,1.05,1.47];
% 
% fig_hdl=figure(2);
% line(AOA_slender_exp,LDratio_slender_exp,'Marker','s','Color','r','LineStyle','none','MarkerFaceColor','r');
% line(AOA_slender_SU2,LDratio_slender_SU2,'Marker','o','Color','b','LineWidth',1,'LineStyle','-');
% line(AOA_list,LDratio_list,'Marker','^','Color',[0.9290 0.6940 0.1250],'LineWidth',1,'LineStyle','--');
% set(gca,'YLim',[-0.2,1.8]);
% xlabel('\alpha/°');
% ylabel('L/D');
% legend('Exp','HF','LF','Location','northwest')
% fig_hdl.set('Position',[488,342,280,210])
% mesh_data on;

% load('blunt_cone_exp.mat');
% load('blunt_cone_SU2.mat')
% load('blunt_cone_HATS.mat')
% 
% 
% fig_hdl=figure();
% line(X_blunt_cone_exp,HF_blunt_cone_exp,'Marker','s','Color','r','LineStyle','none','MarkerFaceColor','r');
% line(X_blunt_cone_SU2/568.7e-3,HF_blunt_cone_SU2/HF_0,'Color','b','LineWidth',1);
% line(X_blunt_cone_HATS/568.7e-3,HF_blunt_cone_HATS/HF_0,'Color',[0.9290 0.6940 0.1250],'LineWidth',1);
% set(gca,'Xlim',[-0.1,1]);
% set(gca,'Ylim',[-0.3,2.3]);
% xlabel('x/L_b');
% ylabel('HF/HF_0');
% legend('Exp','HF','LF','Location','northeast')
% fig_hdl.set('Position',[488,342,420,315])



