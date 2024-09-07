clc;
clear;
close all;

%% PATH

% addpath cfg
% addpath input
% 
% addpath mesh
% addpath src_base
% addpath src_solver
% 
% addpath cgns4m
% startup_cgns4m

%% run

% [model,CA]=PanATS('slender.cfg');
[model,CA]=PanATS('blunt_cone.cfg');
% [model,CA]=PanATS('hermes.cfg'); 
% [model,CA]=PanATS('sanger.cfg');
% [model,CA]=PanATS('HL20_Ma6.cfg');
% [model,CA]=PanATS('HL20_Ma10.cfg');
% [model,CA]=PanATS('waverider.cfg'); 

% config=PanATSConfig('airfoil.cfg');
% config=config.data_dict;

% % cylinder flow
% num=160;
% i=(1:num)-0.5;i=(i/num*2*pi)';
% pnt_list=[1-cos(i),sin(i)];
% elem_list=[[num,1:num-1];1:num]';
% geometry.point_list=pnt_list;
% geometry.dimension=2;
% airfoil.element_list=elem_list;
% airfoil.ID=3;
% airfoil.number_list=2;
% mesh_data.geometry=geometry;
% mesh_data.airfoil=airfoil;
% 
% config.mesh_data=mesh_data;
% [model,CA]=PanATS(config);

%% display

% displayModel([],model,'Cp')
% displayModel([],model,'SL')
% displayModel([],model,'HF')
% displayModel([],model,'Cf')
