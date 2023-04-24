clc;
clear;
close all hidden;

addpath([pwd,'\src']);
addpath([pwd,'\input']);
addpath([pwd,'\mesh']);
addpath([pwd,'\lib']);

global user_model

% user_model=preModelCFG('slender.cfg');
% user_model=preModelCFG('blunt_cone.cfg');
% user_model=preModelCFG('waverider.cfg');
% user_model=preModelCFG('INP_plate.cfg');
% user_model=preModelCFG('geo_test.cfg');

preModelPanel();
[Cl,Cd,LDratio,Cx,Cy,Cz,Cmx,Cmy,Cmz]=solveModelHypersonicInviscid();
[max_streamline_len]=solveModelStreamline();
[max_heat_flow]=solveModelHypersonicHeat();
[Cl,Cd,LDratio,Cx,Cy,Cz,Cmx,Cmy,Cmz]=solveModelHypersonicViscid();

% [area,volume]=solveGeometry()
% displayMarker('SLENDER')
% displayMarker('BLUNT')
% displayMarker('waverider')
% displayMarker('Part-plate')

% displayModel('Cp')
% displayModel('SL')
% displayModel('Q')
% displayModel('Cf')


