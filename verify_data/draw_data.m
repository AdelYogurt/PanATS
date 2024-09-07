clc;
clear;
close all hidden;

%% slender

% load('slender_exp.mat');
% load('slender_SU2.mat')
% load('slender_PanATS.mat')
% 
% fig_hdl_aero=figure();
% 
% % CD
% axe_hdl_CD=axes(fig_hdl_aero);
% line(AOA_exp,CD_exp,'Marker','^','Color','k','MarkerFaceColor','k','LineStyle','None');
% line(AOA_SU2,CD_SU2,'Marker','s','Color','r');
% line(AOA_PanATS,CD_PanATS,'Marker','o','Color','b');
% grid on;box on;
% xlabel('\fontname{times new roman}\alpha/°','FontSize',10.5);
% ylabel('\fontname{times new roman}C_{D}','FontSize',10.5);
% 
% % CL
% axe_hdl_CL=axes(fig_hdl_aero);
% line(AOA_exp,CL_exp,'Marker','^','Color','k','MarkerFaceColor','k','LineStyle','None');
% line(AOA_SU2,CL_SU2,'Marker','s','Color','r');
% line(AOA_PanATS,CL_PanATS,'Marker','o','Color','b');
% grid on;box on;
% xlabel('\fontname{times new roman}\alpha/°','FontSize',10.5);
% ylabel('\fontname{times new roman}C_{L}','FontSize',10.5);
% 
% % CEff
% axe_hdl_CEff=axes(fig_hdl_aero);
% line(AOA_exp,CEff_exp,'Marker','^','Color','k','MarkerFaceColor','k','LineStyle','None');
% line(AOA_SU2,CEff_SU2,'Marker','s','Color','r');
% line(AOA_PanATS,CEff_PanATS,'Marker','o','Color','b');
% grid on;box on;
% xlabel('\fontname{times new roman}\alpha/°','FontSize',10.5);
% ylabel('\fontname{times new roman}C_{L/D}','FontSize',10.5);
% 
% % CMy
% axe_hdl_CMy=axes(fig_hdl_aero);
% exp_li=line(AOA_exp,CMy_exp,'Marker','^','Color','k','MarkerFaceColor','k','LineStyle','None');
% SU2_li=line(AOA_SU2,CMy_SU2,'Marker','s','Color','r');
% PanATS_li=line(AOA_PanATS,CMy_PanATS,'Marker','o','Color','b');
% grid on;box on;
% xlabel('\fontname{times new roman}\alpha/°','FontSize',10.5);
% ylabel('\fontname{times new roman}C_{My}','FontSize',10.5);
% 
% % figure
% axe_hdl_CD.set('Position',[0.15,0.65,0.32,0.24],'FontSize',10.5)
% axe_hdl_CL.set('Position',[0.65,0.65,0.32,0.24],'FontSize',10.5)
% axe_hdl_CEff.set('Position',[0.15,0.20,0.32,0.24],'FontSize',10.5)
% axe_hdl_CMy.set('Position',[0.65,0.20,0.32,0.24],'FontSize',10.5)
% 
% annotation_hdl=annotation('textbox',[0.14,0.52,0.38,0.01],'EdgeColor','none','FontSize',10.5,'String','\fontname{Times New Roman}(a)  \fontname{宋体}阻力系数');
% annotation_hdl=annotation('textbox',[0.64,0.52,0.38,0.01],'EdgeColor','none','FontSize',10.5,'String','\fontname{Times New Roman}(b)  \fontname{宋体}升力系数');
% annotation_hdl=annotation('textbox',[0.12,0.08,0.38,0.01],'EdgeColor','none','FontSize',10.5,'String','\fontname{Times New Roman}(c)  \fontname{宋体}升阻比系数');
% annotation_hdl=annotation('textbox',[0.60,0.08,0.42,0.01],'EdgeColor','none','FontSize',10.5,'String','\fontname{Times New Roman}(d)  \fontname{宋体}俯仰力矩系数');
% 
% lgd_hdl=legend([exp_li,SU2_li,PanATS_li],{'\fontname{宋体}实验数据','\fontname{宋体}SU2','\fontname{宋体}PanATS'});
% lgd_hdl.set('Parent',fig_hdl_aero,'Orientation','horizontal','Position',[0.10,0.92,0.80,0.06],'FontSize',10.5,'box','off');
% fig_hdl_aero.set('Position',[500,350,340,380]);
% 
% % print(fig_hdl_aero,'slender_verify.emf','-dmeta','-r600');
% % print(fig_hdl_aero,'slender_verify.png','-dpng','-r1200');

%% blunt cone

load('blunt_cone_exp.mat');
% load('blunt_cone_SU2.mat');
load('blunt_cone_PanATS.mat');

fig_hdl=figure();axe_hdl=axes(fig_hdl);
line(axe_hdl,X_exp,HF_exp,'Marker','s','Color','r','LineStyle','none','MarkerFaceColor','r');
% line(axe_hdl,X_SU2/568.7e-3,HF_SU2/HF_0,'Color','b','LineWidth',1);
line(axe_hdl,X_PanATS/568.7e-3,HF_PanATS/HF_0,'Color',[0.9290 0.6940 0.1250],'LineWidth',1);
axe_hdl.set('FontSize',12,'FontName','times new roman');
xlim([-0.1,1]);
ylim([-0.1,0.4]);
xlabel('\itX/L_b','FontSize',12,'FontName','times new roman');
ylabel('\itQ/Q_0','FontSize',12,'FontName','times new roman');
legend('Exp','PanATS','Location','northeast','FontSize',12,'FontName','times new roman')
fig_hdl.set('Position',[488,342,420,315])
box on;grid on;

%% hermes

% load('hermes_exp_Ma8_aero.mat');
% load('hermes_PanATS_Ma8_aero.mat');
% load('hermes_SU2_Ma8_aero.mat');
% 
% fig_hdl_aero=figure();
% 
% % CD
% axe_hdl_CD=axes(fig_hdl_aero);
% line(AOA_exp_JB,CD_exp_JB,'Marker','^','Color','k','MarkerFaceColor','k','LineStyle','None');
% line(AOA_SU2,CD_SU2,'Marker','s','Color','r');
% line(AOA_PanATS,CD_PanATS,'Marker','o','Color','b');
% grid on;box on;
% xlabel('\fontname{times new roman}\alpha/°','FontSize',10.5);
% ylabel('\fontname{times new roman}C_{D}','FontSize',10.5);
% 
% % CL
% axe_hdl_CL=axes(fig_hdl_aero);
% line(AOA_exp_JB,CL_exp_JB,'Marker','^','Color','k','MarkerFaceColor','k','LineStyle','None');
% line(AOA_SU2,CL_SU2,'Marker','s','Color','r');
% line(AOA_PanATS,CL_PanATS,'Marker','o','Color','b');
% grid on;box on;
% xlabel('\fontname{times new roman}\alpha/°','FontSize',10.5);
% ylabel('\fontname{times new roman}C_{L}','FontSize',10.5);
% 
% % CEff
% axe_hdl_CEff=axes(fig_hdl_aero);
% line(AOA_exp_JB,CEff_exp_JB,'Marker','^','Color','k','MarkerFaceColor','k','LineStyle','None');
% line(AOA_SU2,CEff_SU2,'Marker','s','Color','r');
% line(AOA_PanATS,CEff_PanATS,'Marker','o','Color','b');
% grid on;box on;
% xlabel('\fontname{times new roman}\alpha/°','FontSize',10.5);
% ylabel('\fontname{times new roman}C_{L/D}','FontSize',10.5);
% 
% % CMy
% axe_hdl_CMy=axes(fig_hdl_aero);
% exp_JB_li=line(AOA_exp_JB,CMy_exp_JB,'Marker','^','Color','k','MarkerFaceColor','k','LineStyle','None');
% SU2_li=line(AOA_SU2,CMy_SU2,'Marker','s','Color','r');
% PanATS_li=line(AOA_PanATS,CMy_PanATS,'Marker','o','Color','b');
% grid on;box on;
% xlabel('\fontname{times new roman}\alpha/°','FontSize',10.5);
% ylabel('\fontname{times new roman}C_{My}','FontSize',10.5);
% 
% % figure
% axe_hdl_CD.set('Position',[0.15,0.65,0.32,0.24],'FontSize',10.5)
% axe_hdl_CL.set('Position',[0.65,0.65,0.32,0.24],'FontSize',10.5)
% axe_hdl_CEff.set('Position',[0.15,0.20,0.32,0.24],'FontSize',10.5)
% axe_hdl_CMy.set('Position',[0.65,0.20,0.32,0.24],'FontSize',10.5)
% 
% annotation_hdl=annotation('textbox',[0.14,0.52,0.38,0.01],'EdgeColor','none','FontSize',10.5,'String','\fontname{Times New Roman}(a)  \fontname{宋体}阻力系数');
% annotation_hdl=annotation('textbox',[0.64,0.52,0.38,0.01],'EdgeColor','none','FontSize',10.5,'String','\fontname{Times New Roman}(b)  \fontname{宋体}升力系数');
% annotation_hdl=annotation('textbox',[0.12,0.08,0.38,0.01],'EdgeColor','none','FontSize',10.5,'String','\fontname{Times New Roman}(c)  \fontname{宋体}升阻比系数');
% annotation_hdl=annotation('textbox',[0.60,0.08,0.42,0.01],'EdgeColor','none','FontSize',10.5,'String','\fontname{Times New Roman}(d)  \fontname{宋体}俯仰力矩系数');
% 
% lgd_hdl=legend([exp_JB_li,SU2_li,PanATS_li],{'\fontname{宋体}实验数据','\fontname{宋体}SU2','\fontname{宋体}PanATS'});
% lgd_hdl.set('Parent',fig_hdl_aero,'Orientation','horizontal','Position',[0.10,0.92,0.80,0.06],'FontSize',10.5,'box','off');
% fig_hdl_aero.set('Position',[500,350,340,380]);
% 
% % print(fig_hdl_aero,'hermes_verify.emf','-dmeta','-r600');
% % print(fig_hdl_aero,'hermes_verify.png','-dpng','-r1200');

%% HL20




