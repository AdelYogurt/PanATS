clc;
clear;
close all hidden;

addpath([pwd,'\src']);
addpath([pwd,'\lib']);
addpath([pwd,'\mesh']);

part=readMeshSTL('waverider_wing.stl');

% displayPart(part);
% writeMeshSTL('test.stl',part);

[part_list,point_list,geometry]=convertSTLToMesh(part,false(1));
displayPart(part_list,point_list);
% writeMeshINP('test.inp',part_list,point_list)

% part_list=readMeshWGS('tmx1242.wgs');

% displayPart(part_list);
% writeMeshWGS('test.wgs',part_list);

% [part_list,point_list,geometry]=convertWGSToMesh(part_list,true(1));
% displayPart(part_list,point_list);
% writeMeshINP('test.inp',part_list,point_list)


% [part_list,point_list,geometry]=readMeshINP('mixed_mesh.inp');
% displayPart(part_list,point_list);
% writeMeshINP('test.inp',part_list,point_list);


% [part_list,point_list,geometry]=readMeshSU2('slender.su2');
% displayPart(part_list,point_list);
% writeMeshINP('test.inp',part_list,point_list);

