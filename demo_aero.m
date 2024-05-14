clc;
clear;
close all hidden;

global user_model

%% PATH

addpath cfg
addpath input
addpath mesh
addpath src_base
addpath src_geo
addpath src_solver

%% pre model

% config=PanATSConfig('hermes.cfg'); 
% AOA_PanATS=[-5,0,5,10,15,20,25,30];

% config=PanATSConfig('slender.cfg'); 
% AOA_PanATS=[0,3,5,7,10];

%% calculate

preModelPanel(config);
user_model.config.INFORMATION=0;

CD_PanATS=zeros(size(AOA_PanATS));
CL_PanATS=zeros(size(AOA_PanATS));
CSF_PanATS=zeros(size(AOA_PanATS));
CFx_PanATS=zeros(size(AOA_PanATS));
CFy_PanATS=zeros(size(AOA_PanATS));
CFz_PanATS=zeros(size(AOA_PanATS));
CMx_PanATS=zeros(size(AOA_PanATS));
CMy_PanATS=zeros(size(AOA_PanATS));
CMz_PanATS=zeros(size(AOA_PanATS));
CEff_PanATS=zeros(size(AOA_PanATS));
times_PanATS=zeros(size(AOA_PanATS));

solveModelHypersonicInviscid();
solveModelStreamline();
solveModelBoundaryLayer();
solveModelHypersonicViscid();
solveModelHypersonicHeat();

for idx=1:length(AOA_PanATS)
    user_model.config.AOA=AOA_PanATS(idx);
    tic;

    [CD,CL,CSF,CFx,CFy,CFz,CMx,CMy,CMz,CEff]=solveModelHypersonicInviscid();
    [max_streamline_len]=solveModelStreamline();
    solveModelBoundaryLayer();
    [CD,CL,CSF,CFx,CFy,CFz,CMx,CMy,CMz,CEff]=solveModelHypersonicViscid();
    [max_heat_flux]=solveModelHypersonicHeat();

    CD_PanATS(idx)=CD;
    CL_PanATS(idx)=CL;
    CSF_PanATS(idx)=CSF;
    CFx_PanATS(idx)=CFx;
    CFy_PanATS(idx)=CFy;
    CFz_PanATS(idx)=CFz;
    CMx_PanATS(idx)=CMx;
    CMy_PanATS(idx)=CMy;
    CMz_PanATS(idx)=CMz;
    CEff_PanATS(idx)=CEff;
    times_PanATS(idx)=toc;
end

flight_condition.MACH_NUMBER= 8.04;
flight_condition.FREESTREAM_TEMPERATURE = 63.89;
flight_condition.FREESTREAM_PRESSURE = 773.48;

% save('PanATS','flight_condition','AOA_PanATS','CD_PanATS','CL_PanATS','CSF_PanATS','CFx_PanATS','CFy_PanATS','CFz_PanATS','CMx_PanATS','CMy_PanATS','CMz_PanATS','CEff_PanATS','times_PanATS')

