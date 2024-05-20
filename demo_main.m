clc;
clear;
close all hidden;

global user_model

%% PATH

% addpath cfg
% addpath input
% 
% addpath mesh
% addpath src_base
% addpath src_geo
% addpath src_solver
% 
% addpath cgns4m
% startup_cgns4m

%% pre model

% config=PanATSConfig('slender.cfg');
% config=PanATSConfig('blunt_cone.cfg');
% config=PanATSConfig('INP_plate.cfg');
% config=PanATSConfig('WWD.cfg');
% config=PanATSConfig('hermes.cfg'); 
% config=PanATSConfig('waverider.cfg'); 
% config=PanATSConfig('sanger.cfg');
% config=PanATSConfig('HL20_Ma10.cfg');
% config=PanATSConfig('sym.cfg');
preModelPanel(config);

%% solve model

[area,area_x,area_y,area_z,volume]=solveGeometry();

% inviscid
[CD,CL,CSF,CFx,CFy,CFz,CMx,CMy,CMz,CEff]=solveModelHypersonicInviscid();

% viscid
[max_streamline_len]=solveModelStreamline();
solveModelBoundaryLayer();
[CD,CL,CSF,CFx,CFy,CFz,CMx,CMy,CMz,CEff]=solveModelHypersonicViscid();
[max_heat_flux]=solveModelHypersonicHeat();

%% post model

postModel()

% displayModel('Cp')
% displayModel('SL')
% displayModel('HF')
% displayModel('Cf')
