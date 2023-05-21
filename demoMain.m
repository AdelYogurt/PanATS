clc;
clear;
% close all hidden;

addpath([pwd,'\src']);
addpath([pwd,'\input']);
addpath([pwd,'\mesh']);
addpath([pwd,'\lib']);

global user_model

% user_model=preModelCFG('slender.cfg');
% user_model=preModelCFG('blunt_cone.cfg');
% user_model=preModelCFG('waverider.cfg');
% user_model=preModelCFG('INP_plate.cfg');
user_model=preModelCFG('WDB.cfg');

preModelPanel();
[Cl,Cd,LDratio,Cx,Cy,Cz,Cmx,Cmy,Cmz]=solveModelHypersonicInviscid();
[max_streamline_len]=solveModelStreamline();
solveModelBoundaryLayer();
[Cl,Cd,LDratio,Cx,Cy,Cz,Cmx,Cmy,Cmz]=solveModelHypersonicViscid();
[max_heat_flux]=solveModelHypersonicHeat();

% [area,volume]=solveGeometry()
% displayMarker('SLENDER')
% displayMarker('BLUNT')
% displayMarker('waverider')
% displayMarker()

% displayModel('Cp')
% displayModel('SL')
displayModel('Q')
% displayModel('Cf')


