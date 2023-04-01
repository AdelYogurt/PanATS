clc;
clear;
close all hidden;

% AOA_list=[0,3,5,7,10];
% Cl_list=zeros(1,length(AOA_list));
% Cd_list=zeros(1,length(AOA_list));
% LDratio_list=zeros(1,length(AOA_list));
% Cmy_list=zeros(1,length(AOA_list));
% for AOA_index=1:length(AOA_list)
%     AOA=AOA_list(AOA_index);
%     g_model.AOA=AOA;
%     [Cl,Cd,LDratio,Cmy]=solveModelHypersonicInviscid();
%     Cl_list(AOA_index)=Cl;
%     Cd_list(AOA_index)=Cd;
%     LDratio_list(AOA_index)=LDratio;
%     Cmy_list(AOA_index)=Cmy;
% end
% 
% figure();
% sgtitle('Slender Aerodynamic Verification');
% 
% subplot(2,2,1);
% line(AOA_list,[0,0.07,0.16,0.31,0.58],'Marker','o','Color','r','LineStyle','none');
% line(AOA_list,Cl_list,'Marker','*','Color','b');
% set(gca,'YLim',[-0.1,0.7]);
% xlabel('\alpha deg');
% ylabel('C_L');
% 
% subplot(2,2,2);
% line(AOA_list,[0.2,0.23,0.24,0.27,0.36],'Marker','o','Color','r','LineStyle','none');
% line(AOA_list,Cd_list,'Marker','*','Color','b');
% set(gca,'YLim',[0.0,0.4]);
% xlabel('\alpha deg');
% ylabel('C_D');
% 
% subplot(2,2,3);
% line(AOA_list,[0,0.33,0.68,1.15,1.6],'Marker','o','Color','r','LineStyle','none');
% line(AOA_list,LDratio_list,'Marker','*','Color','b');
% set(gca,'YLim',[-0.2,1.8]);
% xlabel('\alpha deg');
% ylabel('L/D');
% 
% subplot(2,2,4);
% line(AOA_list,[0,0.03,0.07,0.095,0.158],'Marker','o','Color','r','LineStyle','none');
% line(AOA_list,Cmy_list,'Marker','*','Color','b')
% set(gca,'YLim',[-0.05,0.2]);
% xlabel('\alpha deg');
% ylabel('C_m');

% preModelPanel(AOA,SIDESLIP_ANGLE,SYMMETRY)
% solveModelStreamline...
%     ([],[],[],[],[],[],AOA,SIDESLIP_ANGLE);
% [Cl,Cd,LDratio,Cmz]=solveModelHypersonicInviscid...
%     ([],[],FREESTREAM_TEMPERATURE,FREESTREAM_PRESSURE,MACH_NUMBER,GAMA_VALUE,AOA,SIDESLIP_ANGLE,...
%     [REF_ORIGIN_MOMENT_X,REF_ORIGIN_MOMENT_Y,REF_ORIGIN_MOMENT_Z],REF_LENGTH,REF_AREA,REYNOLDS_NUMBER);
% max_heat_flow=solveModelHypersonicHeat...
%     ([],[],FREESTREAM_TEMPERATURE,FREESTREAM_PRESSURE,MACH_NUMBER,GAMA_VALUE,AOA,SIDESLIP_ANGLE,...
%     MARKER_ISOTHERMAL);
% [Cl,Cd,LDratio,Cmz]=solveModelHypersonicViscid...
%     ([],[],FREESTREAM_TEMPERATURE,FREESTREAM_PRESSURE,MACH_NUMBER,GAMA_VALUE,AOA,SIDESLIP_ANGLE,...
%     [REF_ORIGIN_MOMENT_X,REF_ORIGIN_MOMENT_Y,REF_ORIGIN_MOMENT_Z],REF_LENGTH,REF_AREA,REYNOLDS_NUMBER)

% Elastic_modulus=153e9;
% Poisson_ratio=0.28;
% density=8.98e3;
% 
% thickness=0.01;
% preModelFEM(FREESTREAM_PRESSURE,SYMMETRY,thickness,Elastic_modulus,Poisson_ratio,density)
% solveModelFEM(1)

% [area,volume]=solveGeometry();