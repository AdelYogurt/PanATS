clc;
clear;
close all hidden;

addpath([pwd,'\cfg']);
addpath([pwd,'\input']);
addpath([pwd,'\mesh']);
addpath([pwd,'\src_base']);
addpath([pwd,'\src_geo']);
addpath([pwd,'\src_solver']);

global user_model

config=PanATSConfig('slender.cfg');
% config=PanATSConfig('blunt_cone.cfg');
% config=PanATSConfig('INP_plate.cfg');
% config=PanATSConfig('WWD.cfg');
% config=PanATSConfig('hermes.cfg'); 
% config=PanATSConfig('waverider.cfg'); 
% config=PanATSConfig('sanger.cfg')
% config=PanATSConfig('HL20.cfg');
preModelPanel(config);

[area,area_x,area_y,area_z,volume]=solveGeometry();

% inviscid
[CL,CD,CEff,CFx,CFy,CFz,CMx,CMy,CMz]=solveModelHypersonicInviscid();

% viscid
[max_streamline_len]=solveModelStreamline();
solveModelBoundaryLayer();
[CL,CD,CEff,CFx,CFy,CFz,CMx,CMy,CMz]=solveModelHypersonicViscid();
[max_heat_flux]=solveModelHypersonicHeat();

displayModel('Cp')
% displayModel('SL')
% displayModel('HF')
% displayModel('Cf')
