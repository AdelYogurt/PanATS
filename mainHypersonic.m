clc;
% clear all;
close all hidden;

addpath([pwd,'\src']);

global g_geometry g_Point g_Element g_Marker...
    ADtree_marker_element HATS_element_list...
    streamline_output inviscid_output heat_output viscid_output FEM_output...
    post_output

% Mach number (non-dimensional, based on the free-stream values)
% MACH_NUMBER=15;
% MACH_NUMBER=13.8;
MACH_NUMBER=8.2;
% MACH_NUMBER=10.6;
% Free-stream pressure (101325.0 N/m^2, 2116.216 psf by default)
FREESTREAM_PRESSURE=951.5;
% FREESTREAM_PRESSURE=132;
% Free-stream temperature (288.15 K by default)
FREESTREAM_TEMPERATURE=89.3;
% FREESTREAM_TEMPERATURE=47.3397;
% Angle of attack (degrees)
AOA=0.0;
% Side-slip angle (degrees, only for compressible flows)
SIDESLIP_ANGLE=0.0;

% Reynolds number (non-dimensional, based on the free-stream values) Reynolds length (1 m by default)
REYNOLDS_NUMBER= 9.35E6;

GAMA_VALUE=1.4;
% Marker(s) of the surface where the functional (Cd, Cl, etc.) will be evaluated
% MARKER_MONITORING={'WAVERIDER'};
% Navier-Stokes wall boundary marker(s) (NONE = no marker) constant wall temperature
MARKER_ISOTHERMAL=300;
% Reference origin for moment computation
REF_ORIGIN_MOMENT_X = 0.125;
REF_ORIGIN_MOMENT_Y = 0.00;
REF_ORIGIN_MOMENT_Z = 0;
% Reference length for pitching, rolling, and yawing non-dimensional moment
REF_LENGTH = 0.177;

% Reference area for force coefficients (0 implies automatic calculation)
REF_AREA = pi*1e-4/2;
% REF_AREA = 1;

% Symmetry faca
% SYMMETRY = 'ZOX';
SYMMETRY = [];

% mesh_filename='cubic';
% [g_geometry,g_Point,g_Element,g_Marker]=readMeshDataSTL(mesh_filename);

% mesh_filename='waverider_cone';
% [g_geometry,g_Point,g_Element,g_Marker]=readMeshDataSTL(mesh_filename,0.001);

% mesh_filename='waverider_temp';
% [g_geometry,g_Point,g_Element,g_Marker]=readMeshDataSTL(mesh_filename);

% mesh_filename='BWB';
% [g_geometry,g_Point,g_Element,g_Marker]=readMeshDataSTL(mesh_filename);

% mesh_filename='waverider_panel';
% [g_geometry,g_Point,g_Element,g_Marker]=readMeshDataSU2(mesh_filename,MARKER_MONITORING);

% MARKER_MONITORING={'SLENDER'};
% mesh_filename='slender.su2';
% [g_geometry,g_Point,g_Element,g_Marker]=readMeshDataSU2(mesh_filename,MARKER_MONITORING);

MARKER_MONITORING={'SLENDER'};
mesh_filename='slender_panel.su2';
[g_geometry,g_Point,g_Element,g_Marker]=readMeshDataSU2(mesh_filename,MARKER_MONITORING);

% MARKER_MONITORING={'BLUNT'};
% mesh_filename='blunt_panel.su2';
% [g_geometry,g_Point,g_Element,g_Marker]=readMeshDataSU2(mesh_filename,MARKER_MONITORING);

% AOA_list=[0,3,5,7,10];
% Cl_list=zeros(1,length(AOA_list));
% Cd_list=zeros(1,length(AOA_list));
% LDratio_list=zeros(1,length(AOA_list));
% Cmz_list=zeros(1,length(AOA_list));
% for AOA_index=1:length(AOA_list)
%     AOA=AOA_list(AOA_index);
%     preModelPanel(AOA,SIDESLIP_ANGLE,SYMMETRY)
%     solveModelStreamline...
%         ([],[],[],[],[],[],AOA,SIDESLIP_ANGLE);
%     [Cl,Cd,LDratio,Cmz]=solveModelHypersonicInviscid...
%         ([],[],FREESTREAM_TEMPERATURE,FREESTREAM_PRESSURE,MACH_NUMBER,GAMA_VALUE,AOA,SIDESLIP_ANGLE,...
%         [REF_ORIGIN_MOMENT_X,REF_ORIGIN_MOMENT_Y,REF_ORIGIN_MOMENT_Z],REF_LENGTH,REF_AREA,REYNOLDS_NUMBER);
%     max_heat_flow=solveModelHypersonicHeat...
%         ([],[],FREESTREAM_TEMPERATURE,FREESTREAM_PRESSURE,MACH_NUMBER,GAMA_VALUE,AOA,SIDESLIP_ANGLE,...
%         MARKER_ISOTHERMAL);
%     [Cl,Cd,LDratio,Cmz]=solveModelHypersonicViscid...
%         ([],[],FREESTREAM_TEMPERATURE,FREESTREAM_PRESSURE,MACH_NUMBER,GAMA_VALUE,AOA,SIDESLIP_ANGLE,...
%         [REF_ORIGIN_MOMENT_X,REF_ORIGIN_MOMENT_Y,REF_ORIGIN_MOMENT_Z],REF_LENGTH,REF_AREA,REYNOLDS_NUMBER);
%     Cl_list(AOA_index)=Cl;
%     Cd_list(AOA_index)=Cd;
%     LDratio_list(AOA_index)=LDratio;
%     Cmz_list(AOA_index)=Cmz;
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
% line(AOA_list,Cmz_list,'Marker','*','Color','b')
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

% displayModel('log_P')
% displayMarker('WAVERIDER')

% rmpath([pwd,'\src']);
