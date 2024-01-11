function solveModelFEM(inertial_release)
% solve model
% process border,calculate stiffness matrix and assembly,solve delta
% calculate stress of element
% notice shear strain is engineering shear strain(twice as much as theoretical mechanics)
% Mindlin plate + film
%
% g_Point is coordinate of all node
% g_Element is element' node
% g_Material is meterial list
% g_Border is border displace
% g_Focus is concentrated force
% g_Focus is face force
% g_Volume is volume force
%
% g_Local: {vector_x;vector_y;vector_z;...}
% g_Point_Local: {x1,y1,t1;...}
%
% delta format:
% {u_1,v_1,w_1,theta_x1,theta_y1,theta_z1,...}
% sigma format:
% {stress_xx,stress_yy,stress_zz,stress_xy,stress_yz,stress_zx,stress_von}
% F format:
% {F_x1,F_y1,F_z1,M_x1,M_y1,M_z1,...}
%
% copyright 2022.12 Adel
%
global g_geometry g_Point g_Element g_Material g_Local g_Point_Local...
    g_Border g_Focus g_Face g_Volume g_K g_M g_Delta...
    FEM_output

INFORMATION_FLAG=0;

reduced_integrate_flag=0;
node_number=3;
DOF_node=6;

if INFORMATION_FLAG
    fprintf('solveModelFEM: start solve initialize\n');
end
% step_0,initialize
[point_number,~]=size(g_Point);
% initialize local point coordinate information
% g_K=sparse(point_number*DOF_node,point_number*DOF_node);
% g_M=sparse(point_number*DOF_node,point_number*DOF_node);
% F=sparse(point_number*DOF_node,1);
g_K=zeros(point_number*DOF_node,point_number*DOF_node);
g_M=zeros(point_number*DOF_node,point_number*DOF_node);
F=zeros(point_number*DOF_node,1);

% step_1,calculate stiffness matrix and assembly
[element_number,~]=size(g_Element);
% initialize local coordinate information
g_Local=zeros(element_number*3,3);
g_Point_Local=zeros(element_number*node_number,2);
for element_index=1:element_number
    stiffness=stiffness3DShellDKTTri3N(element_index,reduced_integrate_flag);
    mass=mass3DShellDKTTri3N(element_index,reduced_integrate_flag);
    assemble3DShellDKTTri3N(element_index,stiffness,mass);
end

% filename_input='Job-1_Matrix_STIF2.mtx';
% stiffness_ABAQUS=readABAQUSStiffness(filename_input);

% step_2,process concentrate force
[focus_number,~]=size(g_Focus);
for focus_index=1:focus_number
    point_index=g_Focus(focus_index,1);
    direct=g_Focus(focus_index,2);
    focus=g_Focus(focus_index,3);
    F_index=(point_index-1)*DOF_node+direct;
    F(F_index)=focus;
end

% step_3,process face force
[face_number,~]=size(g_Face);
for face_index=1:face_number
    equivalent_force=face3DShellDKTTri3N...
        (g_Face(face_index,:));
    node_index=g_Element(g_Face(face_index,1),2:4);
    F(DOF_node*(node_index-1)+1)=F(DOF_node*(node_index-1)+1)+equivalent_force(:,1);
    F(DOF_node*(node_index-1)+2)=F(DOF_node*(node_index-1)+2)+equivalent_force(:,2);
    F(DOF_node*(node_index-1)+3)=F(DOF_node*(node_index-1)+3)+equivalent_force(:,3);
    F(DOF_node*(node_index-1)+4)=F(DOF_node*(node_index-1)+4)+equivalent_force(:,4);
    F(DOF_node*(node_index-1)+5)=F(DOF_node*(node_index-1)+5)+equivalent_force(:,5);
    F(DOF_node*(node_index-1)+6)=F(DOF_node*(node_index-1)+6)+equivalent_force(:,6);
end

% step_4,process volume force
[volume_number,~]=size(g_Volume);
for volume_index=1:volume_number
    equivalent_force=volume3DShellDKTTri3N...
        (g_Volume(volume_index,:));
    node_index=g_Element(g_Volume(volume_index,1),2:4);
    F(DOF_node*(node_index-1)+1)=F(DOF_node*(node_index-1)+1)+equivalent_force(:,1);
    F(DOF_node*(node_index-1)+2)=F(DOF_node*(node_index-1)+2)+equivalent_force(:,2);
    F(DOF_node*(node_index-1)+3)=F(DOF_node*(node_index-1)+3)+equivalent_force(:,3);
    F(DOF_node*(node_index-1)+4)=F(DOF_node*(node_index-1)+4)+equivalent_force(:,4);
    F(DOF_node*(node_index-1)+5)=F(DOF_node*(node_index-1)+5)+equivalent_force(:,5);
    F(DOF_node*(node_index-1)+6)=F(DOF_node*(node_index-1)+6)+equivalent_force(:,6);
end

% step 4s, process Inertial release
if inertial_release
    [mode_vector,frequency]=eig(g_M\g_K);
    frequency=sqrt(frequency);
    mode_vector_lamda=(mode_vector'*g_M)\(mode_vector'*F);
    F_inr=real(g_M*mode_vector_lamda);
    F=F-F_inr;
end

% step_5,process border displacement,set one method
g_K_solve=g_K;
[border_number,~]=size(g_Border);
for border_index=1:border_number
    point_index=g_Border(border_index,1);
    direct=g_Border(border_index,2);
    displace=g_Border(border_index,3);
    F_index=(point_index-1)*DOF_node+direct;
    F([1:F_index-1,F_index+1:end])=F([1:F_index-1,F_index+1:end])-...
        displace*g_K_solve([1:F_index-1,F_index+1:end],F_index);
    F(F_index)=displace;
    g_K_solve(:,F_index)=0;
    g_K_solve(F_index,:)=0;
    g_K_solve(F_index,F_index)=1;
end
% % step_5,process border displacement,times large number method
% g_K_solve=g_K;
% [border_number,~]=size(g_Border);
% for border_index=1:1:border_number
%     point_index=g_Border(border_index,1);
%     direct=g_Border(border_index,2);
%     displace=g_Border(border_index,3);
%     F_index=(point_index-1)*DOF_node+direct;
%     F(F_index)=displace*g_K_solve(F_index,F_index)*1e12;
%     g_K_solve(F_index,F_index)=g_K_solve(F_index,F_index)*1e12;
% end

% step_6,solve equation
if INFORMATION_FLAG
    fprintf('solveModelFEM: start solve model\n');
end
warning off;
g_Delta=mldivide(g_K_solve,F);
g_Delta=(g_Delta);
warning on;
R=(g_K)*g_Delta-F;

% step_7,calculate surface stress

% element surface stess,format:
% point_index_of_element,% {stress_xx,stress_yy,stress_zz,stress_xy,stress_yz,stress_zx,stress_von},element_index
g_SurfStress=zeros(node_number,7,element_number);

% point surface stress,format:
% % {stress_xx,stress_yy,stress_zz,stress_xy,stress_yz,stress_zx,stress_von},repeat_times
g_SurfStress_Point=zeros(point_number,8);

for element_index=1:element_number
    node_index=g_Element(element_index,2:1+node_number);
    g_SurfStress(:,:,element_index)=stress3DShellDKTTri3N(element_index,reduced_integrate_flag);
    % add point stress
    g_SurfStress_Point(node_index,1:7)=...
        g_SurfStress_Point(node_index,1:7)+g_SurfStress(:,:,element_index);
    g_SurfStress_Point(node_index,8)=...
        g_SurfStress_Point(node_index,8)+1;
end
g_SurfStress_Point(:,1:7)=g_SurfStress_Point(:,1:7)./g_SurfStress_Point(:,8);

if INFORMATION_FLAG
    fprintf('solveModelFEM: successfully solve model\n');
end
    
FEM_output.U_list=g_Delta;
FEM_output.R_list=R;
FEM_output.surf_stress_list=g_SurfStress;
FEM_output.surf_stress_point_list=g_SurfStress_Point;
end

function stiffness=stiffness3DShellDKTTri3N...
    (element_index,reduced_integrate_flag)
% calculate the stiffness matrix of the element
% for thick shell,use gauss integral
% thick direction use quintic Simpson or cubic gaussian integral
%
% delta format:
% {u_1,v_1,w_1,theta_x1,theta_y1,...}
%
global g_Point g_Element g_Local g_Point_Local

% one Gaussian integral point and weight
x1_list=[0.3333333333333333,0.3333333333333333,0.3333333333333333];
w1_list=1;
% quadratic Gaussian integral point and weight
% x2_list=[
%     0.3333333333333333,0.3333333333333333,0.3333333333333333;
%     0.0000000000000000,0.5000000000000000,0.5000000000000000;
%     0.5000000000000000,0.5000000000000000,0.0000000000000000;
%     0.5000000000000000,0.0000000000000000,0.5000000000000000;];
% w2_list=[0.5;0.1666666666666667;0.1666666666666667;0.1666666666666667];

% x2_list=[
%     0.0000000000000000,0.5000000000000000,0.5000000000000000;
%     0.5000000000000000,0.5000000000000000,0.0000000000000000;
%     0.5000000000000000,0.0000000000000000,0.5000000000000000;];
% w2_list=[0.3333333333333333;0.3333333333333333;0.3333333333333333];

x2_list=[
    0.6666666666666667,0.1666666666666667,0.1666666666666667;
    0.1666666666666667,0.6666666666666667,0.1666666666666667;
    0.1666666666666667,0.1666666666666667,0.6666666666666667;];
w2_list=[0.3333333333333333;0.3333333333333333;0.3333333333333333];

% quintic Simpson
x5s_list=[-1,              -0.5,            0,               0.5,             1];
w5s_list=[0.166666666666667,0.666666666666667,0.333333333333333,0.666666666666667,0.166666666666667];

node_number=3;
DOF_node=6;
DOF_element=DOF_node*node_number;

node_index=g_Element(element_index,2:4);
x_list=g_Point(node_index,1);
y_list=g_Point(node_index,2);
z_list=g_Point(node_index,3);
t_list=g_Point(node_index,4);

% calculate shell element geometry parameter
% calculate local coordinate z(normal vector)
e12_x=x_list(2)-x_list(1);e12_y=y_list(2)-y_list(1);e12_z=z_list(2)-z_list(1);
e23_x=x_list(3)-x_list(2);e23_y=y_list(3)-y_list(2);e23_z=z_list(3)-z_list(2);
vector_z=[
    e12_y*e23_z-e23_y*e12_z,...
    e23_x*e12_z-e12_x*e23_z,...
    e12_x*e23_y-e12_y*e23_x];
vector_z=vector_z/norm(vector_z);

% calculate local coordinate y(point from point 1 to point 2)
center_point=[sum(x_list)/4,sum(y_list)/4,sum(z_list)/4];
vector_y=[x_list(2)-x_list(1),y_list(2)-y_list(1),z_list(2)-z_list(1)];
vector_y=vector_y/sqrt(sum(vector_y.^2));

% calculate local coordinate x(vector_y cross vector_z)
vector_x=[
    vector_y(2)*vector_z(3)-vector_y(3)*vector_z(2),...
    vector_y(3)*vector_z(1)-vector_y(1)*vector_z(3),...
    vector_y(1)*vector_z(2)-vector_y(2)*vector_z(1)];

g_Local((element_index*3-2):element_index*3,:)=[vector_x;vector_y;vector_z];

% rank_vector*[vector_x;vector_y;vector_z]' means global to local
% [vector_x;vector_y;vector_z]*colume_vector means global to local

% project x,y,z to local coordinate to get x_local,y_local
coordinate_local=([x_list,y_list,z_list]-center_point)*[vector_x;vector_y;vector_z]';
x_local=coordinate_local(:,1);
y_local=coordinate_local(:,2);
t_local=t_list;

g_Point_Local(node_number*element_index-2:node_number*element_index,1)=x_local;
g_Point_Local(node_number*element_index-2:node_number*element_index,2)=y_local;
g_Point_Local(node_number*element_index-2:node_number*element_index,3)=t_local;

D=calMatrixD(element_index);
stiffness_local=zeros(DOF_element);
stiffness=zeros(DOF_element);

% integral in local coordinate so need to times A and dz_dzeta
if reduced_integrate_flag
    for k=1:length(x5s_list) % zeta
        for integral_index=1:size(x1_list,1)
            x=x1_list(integral_index,:)*x_local;
            y=x1_list(integral_index,:)*y_local;
            % calculate shearing stress,reduce intergrate
            [B,A,dz_dzeta]=calMatrixB...
                (element_index,x,y,x5s_list(k));
            stiffness_local=stiffness_local+...
                w1_list(integral_index)*w5s_list(k)*B'*D*B*dz_dzeta*A;
        end
    end
else
    for k=1:length(x5s_list) % zeta
        for integral_index=1:size(x2_list,1)
            x=x2_list(integral_index,:)*x_local;
            y=x2_list(integral_index,:)*y_local;
            % calculate bending stress
            [B,A,dz_dzeta]=calMatrixB...
                (element_index,x,y,x5s_list(k));
            stiffness_local=stiffness_local+...
                w2_list(integral_index)*w5s_list(k)*B'*D*B*dz_dzeta*A;
        end
    end
end

lamada=[[vector_x;vector_y;vector_z],zeros(3);zeros(3),[vector_x;vector_y;vector_z]];
% project local to global
for rank_index=1:node_number
    for colume_index=1:node_number
        stiffness((rank_index-1)*DOF_node+1:rank_index*DOF_node,(colume_index-1)*DOF_node+1:colume_index*DOF_node)=...
            lamada'*stiffness_local((rank_index-1)*DOF_node+1:rank_index*DOF_node,(colume_index-1)*DOF_node+1:colume_index*DOF_node)*lamada;
    end
end
end
function mass=mass3DShellDKTTri3N...
    (element_index,reduced_integrate_flag)
% calculate the stiffness matrix of the element
% for thick shell,use gauss integral
% thick direction use quintic Simpson or cubic gaussian integral
%
% delta format:
% {u_1,v_1,w_1,theta_x1,theta_y1,...}
%
global g_Point g_Element g_Local g_Point_Local g_Material
rou=g_Material(g_Element(element_index,5),3);
% one Gaussian integral point and weight
x1_list=[0.3333333333333333,0.3333333333333333,0.3333333333333333];
w1_list=1;
% quadratic Gaussian integral point and weight
% x2_list=[
%     0.3333333333333333,0.3333333333333333,0.3333333333333333;
%     0.0000000000000000,0.5000000000000000,0.5000000000000000;
%     0.5000000000000000,0.5000000000000000,0.0000000000000000;
%     0.5000000000000000,0.0000000000000000,0.5000000000000000;];
% w2_list=[0.5;0.1666666666666667;0.1666666666666667;0.1666666666666667];

% x2_list=[
%     0.0000000000000000,0.5000000000000000,0.5000000000000000;
%     0.5000000000000000,0.5000000000000000,0.0000000000000000;
%     0.5000000000000000,0.0000000000000000,0.5000000000000000;];
% w2_list=[0.3333333333333333;0.3333333333333333;0.3333333333333333];

x2_list=[
    0.6666666666666667,0.1666666666666667,0.1666666666666667;
    0.1666666666666667,0.6666666666666667,0.1666666666666667;
    0.1666666666666667,0.1666666666666667,0.6666666666666667;];
w2_list=[0.3333333333333333;0.3333333333333333;0.3333333333333333];

% quintic Simpson
x5s_list=[-1,              -0.5,            0,               0.5,             1];
w5s_list=[0.166666666666667,0.666666666666667,0.333333333333333,0.666666666666667,0.166666666666667];

node_number=3;
DOF_node=6;
DOF_element=DOF_node*node_number;

node_index=g_Element(element_index,2:4);
x_list=g_Point(node_index,1);
y_list=g_Point(node_index,2);
z_list=g_Point(node_index,3);
t_list=g_Point(node_index,4);

% calculate shell element geometry parameter
% calculate local coordinate z(normal vector)
e12_x=x_list(2)-x_list(1);e12_y=y_list(2)-y_list(1);e12_z=z_list(2)-z_list(1);
e23_x=x_list(3)-x_list(2);e23_y=y_list(3)-y_list(2);e23_z=z_list(3)-z_list(2);
vector_z=[
    e12_y*e23_z-e23_y*e12_z,...
    e23_x*e12_z-e12_x*e23_z,...
    e12_x*e23_y-e12_y*e23_x];
vector_z=vector_z/norm(vector_z);

% calculate local coordinate y(point from point 1 to point 2)
center_point=[sum(x_list)/4,sum(y_list)/4,sum(z_list)/4];
vector_y=[x_list(2)-x_list(1),y_list(2)-y_list(1),z_list(2)-z_list(1)];
vector_y=vector_y/sqrt(sum(vector_y.^2));

% calculate local coordinate x(vector_y cross vector_z)
vector_x=[
    vector_y(2)*vector_z(3)-vector_y(3)*vector_z(2),...
    vector_y(3)*vector_z(1)-vector_y(1)*vector_z(3),...
    vector_y(1)*vector_z(2)-vector_y(2)*vector_z(1)];

g_Local((element_index*3-2):element_index*3,:)=[vector_x;vector_y;vector_z];

% rank_vector*[vector_x;vector_y;vector_z]' means global to local
% [vector_x;vector_y;vector_z]*colume_vector means global to local

% project x,y,z to local coordinate to get x_local,y_local
coordinate_local=([x_list,y_list,z_list]-center_point)*[vector_x;vector_y;vector_z]';
x_local=coordinate_local(:,1);
y_local=coordinate_local(:,2);
t_local=t_list;

g_Point_Local(node_number*element_index-2:node_number*element_index,1)=x_local;
g_Point_Local(node_number*element_index-2:node_number*element_index,2)=y_local;
g_Point_Local(node_number*element_index-2:node_number*element_index,3)=t_local;

D=calMatrixD(element_index);
mass_local=zeros(DOF_element);
mass=zeros(DOF_element);

% integral in local coordinate so need to times A and dz_dzeta
if reduced_integrate_flag
    for k=1:length(x5s_list) % zeta
        for integral_index=1:size(x1_list,1)
            x=x1_list(integral_index,:)*x_local;
            y=x1_list(integral_index,:)*y_local;
            % calculate shearing stress,reduce intergrate
            [N,A,dz_dzeta]=calMatrixB...
                (element_index,x,y,x5s_list(k));
            mass_local=mass_local+...
                w1_list(integral_index)*w5s_list(k)*N'*rou*N*dz_dzeta*A;
        end
    end
else
    for k=1:length(x5s_list) % zeta
        for integral_index=1:size(x2_list,1)
            x=x2_list(integral_index,:)*x_local;
            y=x2_list(integral_index,:)*y_local;
            % calculate bending stress
            [N,A,dz_dzeta]=calMatrixN...
                (element_index,x,y,x5s_list(k));
            mass_local=mass_local+...
                w2_list(integral_index)*w5s_list(k)*N'*rou*N*dz_dzeta*A;
        end
    end
end

lamada=[[vector_x;vector_y;vector_z],zeros(3);zeros(3),[vector_x;vector_y;vector_z]];
% project local to global
for rank_index=1:node_number
    for colume_index=1:node_number
        mass((rank_index-1)*DOF_node+1:rank_index*DOF_node,(colume_index-1)*DOF_node+1:colume_index*DOF_node)=...
            lamada'*mass_local((rank_index-1)*DOF_node+1:rank_index*DOF_node,(colume_index-1)*DOF_node+1:colume_index*DOF_node)*lamada;
    end
end
end
function assemble3DShellDKTTri3N(element_index,stiffness,mass)
% integrate the element stiffness matrix into the overall stiffness matrix
%
global g_Element g_K g_M
node_number=3;
DOF_node=6;
node_index=repmat(g_Element(element_index,2:1+node_number),DOF_node,1)*DOF_node;
node_index(1,:)=node_index(1,:)-5;
node_index(2,:)=node_index(2,:)-4;
node_index(3,:)=node_index(3,:)-3;
node_index(4,:)=node_index(4,:)-2;
node_index(5,:)=node_index(5,:)-1;

g_K(node_index,node_index)=g_K(node_index,node_index)+stiffness;
g_M(node_index,node_index)=g_M(node_index,node_index)+mass;
end

function [N,A,dz_dzeta]=calMatrixN...
    (element_index,x,y,zeta)
% calculate the strain matrix B of this element
%
global g_Point_Local
node_number=3;
DOF_node=6;
DOF_element=node_number*DOF_node;

x_local=g_Point_Local(node_number*(element_index-1)+1:node_number*element_index,1);
y_local=g_Point_Local(node_number*(element_index-1)+1:node_number*element_index,2);
t_local=g_Point_Local(node_number*(element_index-1)+1:node_number*element_index,3);

% calculate partial derivative of shape function
A=det([[1;1;1],x_local,y_local])/2;

a=[x_local(2)*y_local(3)-x_local(3)*y_local(2),...
    x_local(3)*y_local(1)-x_local(1)*y_local(3),...
    x_local(1)*y_local(2)-x_local(2)*y_local(1)];
a_A=a/A/2;
b1=y_local(2)-y_local(3);
b2=y_local(3)-y_local(1);
b3=y_local(1)-y_local(2);
b_A=[b1,b2,b3]/A/2;
c1=x_local(3)-x_local(2);
c2=x_local(1)-x_local(3);
c3=x_local(2)-x_local(1);
c_A=[c1,c2,c3]/A/2;

% N
L=a_A+b_A*x+c_A*y;
L1=L(1);
L2=L(2);
L3=L(3);

dz_dzeta=0.5*L*t_local;
z=L*t_local*zeta*0.5;

Nb=[
    L1+L1*L1*L2+L1*L1*L3-L1*L2*L2-L1*L3*L3,...
    (b2-b3)*L1*L1*L3+0.5*(b2-b3)*L1*L2*L3,...
    (c2-c3)*L1*L1*L3+0.5*(c2-c3)*L1*L2*L3,...
    ...
    L2+L2*L2*L3+L2*L2*L1-L2*L3*L3-L2*L1*L1,...
    (b3-b1)*L2*L2*L1+0.5*(b3-b1)*L2*L3*L1,...
    (c3-c1)*L2*L2*L1+0.5*(c3-c1)*L2*L3*L1,...
    ...
    L3+L3*L3*L1+L3*L3*L2-L3*L1*L1-L3*L2*L2,...
    (b1-b2)*L3*L3*L2+0.5*(b1-b2)*L3*L1*L2,...
    (c1-c2)*L3*L3*L2+0.5*(c1-c2)*L3*L1*L2];

% calculating the derivative of shape function to global coordinates
% Nu,Nv
Nu=[
    L(1)*(b3*L(2)-b2*L(3))/2,...
    L(2)*(b1*L(3)-b3*L(1))/2,...
    L(3)*(b2*L(1)-b1*L(2))/2];
Nv=[
    L(1)*(c3*L(2)-c2*L(3))/2,...
    L(2)*(c1*L(3)-c3*L(1))/2,...
    L(3)*(c2*L(1)-c1*L(2))/2];

N=zeros(node_number,DOF_element);
% strain matrix
for node_index=1:node_number
    N(:,(DOF_node*(node_index-1)+1):DOF_node*node_index)=[
        L(node_index)	0               0                   0                   0                   Nu(node_index);
        0               L(node_index)   0                   0                   0                   Nv(node_index);
        0               0               Nb(3*node_index-2)	Nb(3*node_index-1)	Nb(3*node_index)	0];
end
end
function [B,A,dz_dzeta]=calMatrixB...
    (element_index,x,y,zeta)
% calculate the strain matrix B of this element
% local x and local y is rotation to local coordinate by normal fo element
%
global g_Element g_Point_Local g_Material
node_number=3;
DOF_node=6;
DOF_element=node_number*DOF_node;

x_local=g_Point_Local(node_number*(element_index-1)+1:node_number*element_index,1);
y_local=g_Point_Local(node_number*(element_index-1)+1:node_number*element_index,2);
t_local=g_Point_Local(node_number*(element_index-1)+1:node_number*element_index,3);

% calculate partial derivative of shape function
A=det([[1;1;1],x_local,y_local])/2;

a=[x_local(2)*y_local(3)-x_local(3)*y_local(2),...
    x_local(3)*y_local(1)-x_local(1)*y_local(3),...
    x_local(1)*y_local(2)-x_local(2)*y_local(1)];
a_A=a/A/2;
b1=y_local(2)-y_local(3);
b2=y_local(3)-y_local(1);
b3=y_local(1)-y_local(2);
b_A=[b1,b2,b3]/A/2;
c1=x_local(3)-x_local(2);
c2=x_local(1)-x_local(3);
c3=x_local(2)-x_local(1);
c_A=[c1,c2,c3]/A/2;

% L
L=a_A+b_A*x+c_A*y;
L1=L(1);
L2=L(2);

dz_dzeta=0.5*L*t_local;
z=L*t_local*zeta*0.5;

% L derivative
dL_dx=b_A;
dL_dy=c_A;

delta=sum(a);
T=[
    b1^2	b2^2	2*b1*b2
    c1^2	c2^2 	2*c1*c2
    2*b1*c1	2*b2*c2	2*(b1*c2+b2*c1) ];

Nii=[
    -4*L2-12*L1+6, ...
    -6*b2*L1+2*b2-3*b2*L2-b3*L2, ...
    2*c2-6*c2*L1-3*c2*L2-c3*L2, ...
    -4*L2, ...
    (-b3+b1)*L2, ...
    -(c3-c1)*L2, ...
    -6+12*L1+8*L2, ...
    b1*L2-6*b2*L1+4*b2-3*b2*L2, ...
    c1*L2+4*c2-6*c2*L1-3*c2*L2 ];
Njj=[
    -4*L1, ...
    -(b2-b3)*L1, ...
    (-c2+c3)*L1, ...
    6-4*L1-12*L2, ...
    b3*L1+6*b1*L2-2*b1+3*b1*L1, ...
    c3*L1+6*c1*L2-2*c1+3*c1*L1, ...
    8*L1-6+12*L2, ...
    6*b1*L2-4*b1+3*b1*L1-b2*L1, ...
    6*c1*L2-4*c1+3*c1*L1-c2*L1];
Nij=[
    -4*L1-4*L2+2, ...
    -3*b2*L1-b3*L1+1/2*b2-b2*L2-1/2*b3+b3*L2, ...
    -3*c2*L1-c3*L1+1/2*c2-c2*L2-1/2*c3+c3*L2, ...
    -4*L1-4*L2+2, ...
    b3*L2+3*b1*L2+1/2*b3-b3*L1-1/2*b1+b1*L1, ...
    c3*L2+3*c1*L2+1/2*c3-c3*L1-1/2*c1+c1*L1, ...
    -4+8*L1+8*L2, ...
    3*b1*L2-3/2*b1+b1*L1-3*b2*L1+3/2*b2-b2*L2, ...
    3*c1*L2-3/2*c1+c1*L1-3*c2*L1+3/2*c2-c2*L2];

Bb=-1/delta^2*T*[Nii; Njj; Nij]*z;

% calculating the derivative of shape function to global coordinates
% Nu,Nv
% Nu=[
%     L(1)*(b3*L(2)-b2*L(3))/2,...
%     L(2)*(b1*L(3)-b3*L(1))/2,...
%     L(3)*(b2*L(1)-b1*L(2))/2];
% Nv=[
%     L(1)*(c3*L(2)-c2*L(3))/2,...
%     L(2)*(c1*L(3)-c3*L(1))/2,...
%     L(3)*(c2*L(1)-c1*L(2))/2];

dNu_dx=[
    dL_dx(1)*(b3*L(2)-b2*L(3))+L(1)*(b3*dL_dx(2)-b2*dL_dx(3)),...
    dL_dx(2)*(b1*L(3)-b3*L(1))+L(2)*(b1*dL_dx(3)-b3*dL_dx(1)),...
    dL_dx(3)*(b2*L(1)-b1*L(2))+L(3)*(b2*dL_dx(1)-b1*dL_dx(2))]/2;
dNu_dy=[
    dL_dy(1)*(b3*L(2)-b2*L(3))+L(1)*(b3*dL_dy(2)-b2*dL_dy(3)),...
    dL_dy(2)*(b1*L(3)-b3*L(1))+L(2)*(b1*dL_dy(3)-b3*dL_dy(1)),...
    dL_dy(3)*(b2*L(1)-b1*L(2))+L(3)*(b2*dL_dy(1)-b1*dL_dy(2))]/2;

dNv_dx=[
    dL_dx(1)*(c3*L(2)-c2*L(3))+L(1)*(c3*dL_dx(2)-c2*dL_dx(3)),...
    dL_dx(2)*(c1*L(3)-c3*L(1))+L(2)*(c1*dL_dx(3)-c3*dL_dx(1)),...
    dL_dx(3)*(c2*L(1)-c1*L(2))+L(3)*(c2*dL_dx(1)-c1*dL_dx(2))]/2;
dNv_dy=[
    dL_dy(1)*(c3*L(2)-c2*L(3))+L(1)*(c3*dL_dy(2)-c2*dL_dy(3)),...
    dL_dy(2)*(c1*L(3)-c3*L(1))+L(2)*(c1*dL_dy(3)-c3*dL_dy(1)),...
    dL_dy(3)*(c2*L(1)-c1*L(2))+L(3)*(c2*dL_dy(1)-c1*dL_dy(2))]/2;

B=zeros(3,DOF_element);

% strain matrix
for node_index=1:node_number
    B(:,(DOF_node*(node_index-1)+1):DOF_node*node_index)=[
        dL_dx(node_index)	0                   Bb(1,(3*(node_index-1)+1):3*node_index)	(dNu_dx(node_index));
        0                   dL_dy(node_index)	Bb(2,(3*(node_index-1)+1):3*node_index)	(dNv_dy(node_index));
        dL_dy(node_index)	dL_dx(node_index)	Bb(3,(3*(node_index-1)+1):3*node_index)	(dNu_dy(node_index)+dNv_dx(node_index));];
end

end
function D=calMatrixD(element_index)
% calculate the elastic matrix D of this element
% use for thick shell
% k equal 6/5 to correct nolinear of stree in thick direction
%
global g_Element g_Material
E=g_Material(g_Element(element_index,5),1);
mu=g_Material(g_Element(element_index,5),2);

A1=mu;
A2=(1-mu)/2;
A3=E/(1-mu^2);

D=A3*[
    1 A1 0;
    A1 1 0;
    0 0 A2;];
end

function equivalent_force=face3DShellDKTTri3N(parameter)
% calculate the equivalent nodal force of the distributed pressure
% notice that only normal pressure is available
% use gaussian integal
%
% return global coordinate equivalent_force(node_number*DOF_node matrix)
%
% equivalent_force format:
% {F_x1,F_y1,F_z1,M_x1,M_y1,M_y1;...}
%
% quadratic Gaussian integral point and weight
global g_Element g_Local g_Point_Local

% one Gaussian integral point and weight
x1_list=[0.3333333333333333,0.3333333333333333,0.3333333333333333];
w1_list=1;
% quadratic Gaussian integral point and weight
% x2_list=[
%     0.3333333333333333,0.3333333333333333,0.3333333333333333;
%     0.0000000000000000,0.5000000000000000,0.5000000000000000;
%     0.5000000000000000,0.5000000000000000,0.0000000000000000;
%     0.5000000000000000,0.0000000000000000,0.5000000000000000;];
% w2_list=[0.5;0.1666666666666667;0.1666666666666667;0.1666666666666667];

x2_list=[
    0.6666666666666667,0.1666666666666667,0.1666666666666667;
    0.1666666666666667,0.6666666666666667,0.1666666666666667;
    0.1666666666666667,0.1666666666666667,0.6666666666666667;];
w2_list=[0.3333333333333333;0.3333333333333333;0.3333333333333333];

% quintic Simpson
x5s_list=[-1,              -0.5,            0,               0.5,             1];
w5s_list=[0.166666666666667,0.666666666666667,0.333333333333333,0.666666666666667,0.166666666666667];

node_number=3;
DOF_node=6;
DOF_element=DOF_node*node_number;

element_index=parameter(1);
face_index=parameter(2);

% obtain face normal vector
vector_x=g_Local(element_index*3-2,:);
vector_y=g_Local(element_index*3-1,:);
vector_z=g_Local(element_index*3,:);

% rank_vector*lamada' means global to local
% lamada*colume_vector means global to local
lamada=[vector_x;vector_y;vector_z];

x_local=g_Point_Local(node_number*(element_index-1)+1:node_number*element_index,1);
y_local=g_Point_Local(node_number*(element_index-1)+1:node_number*element_index,2);
t_local=g_Point_Local(node_number*(element_index-1)+1:node_number*element_index,3);

b=[y_local(2)-y_local(3),y_local(3)-y_local(1),y_local(1)-y_local(2)];
c=[x_local(3)-x_local(2),x_local(1)-x_local(3),x_local(2)-x_local(1)];
switch length(parameter)
    case 3
        % pressure load
        pressure=parameter(3);
        
        switch face_index
            case 1
                pressure_local=pressure;
                face_node_index=[1,2,3];
                EDGE_FLAG=0;
            case 2
                pressure_local=-pressure;
                face_node_index=[1,2,3];
                EDGE_FLAG=0;
            case 3
                bc_index=3;
                face_node_index=[1,2];
                EDGE_FLAG=1;
            case 4
                bc_index=1;
                face_node_index=[2,3];
                EDGE_FLAG=2;
            case 5
                bc_index=2;
                face_node_index=[3,1];
                EDGE_FLAG=1;
            otherwise
                error('face3DShellDKTTri3N: error face index');
        end
        
        if EDGE_FLAG
            % edge load
            delta_x=x_local(face_node_index(2))-x_local(face_node_index(1));
            delta_y=y_local(face_node_index(2))-y_local(face_node_index(1));
            len_edge=norm([delta_x,delta_y]);
            
            pressure_local=[-delta_y,delta_x]*pressure; % colume vector
            
            equivalent_force=zeros(node_number,DOF_number);
            equivalent_force(face_node_index(1),[1,2,6])=...
                (pressure_local*...
                [len_edge/2,0;
                0,len_edge/2;
                b(bc_index)/12*len_edge,c(bc_index)/12*len_edge]'...
                )*t_local(face_node_index(1));
            equivalent_force(face_node_index(2),[1,2,6])=...
                (pressure_local*...
                [len_edge/2,0;
                0,len_edge/2;
                -b(bc_index)/12*len_edge,-c(bc_index)/12*len_edge]'...
                )*t_local(face_node_index(1));
            
            % project to global coordinate
            equivalent_force(face_node_index(1),:)=equivalent_force(face_node_index(1),:)*[lamada,zeros(3);zeros(3),lamada];
            equivalent_force(face_node_index(2),:)=equivalent_force(face_node_index(2),:)*[lamada,zeros(3);zeros(3),lamada];
        else
            % surface load
            A=det([[1;1;1],x_local,y_local])/2;
            
            equivalent_force=pressure_local*A*[
                0,0,1/3,(b(2)-b(3))/24,(c(2)-c(3))/24,0;
                0,0,1/3,(b(3)-b(1))/24,(c(3)-c(1))/24,0;
                0,0,1/3,(b(1)-b(2))/24,(c(1)-c(2))/24,0];
            
            % project to global coordinate
            equivalent_force=[equivalent_force(:,1:3)*lamada,equivalent_force(:,4:6)*lamada];
        end
    case 5
        % surface stress load
        vector_pressure=parameter(3:end);
        switch face_index
            case 1
                zeta=-1;
                EDGE_FLAG=0;
            case 2
                zeta=1;
                EDGE_FLAG=0;
            case 3
                face_node_index=[1,2];
                bc_index=3;
                EDGE_FLAG=1;
            case 4
                face_node_index=[2,3];
                bc_index=1;
                EDGE_FLAG=1;
            case 5
                face_node_index=[3,1];
                bc_index=2;
                EDGE_FLAG=1;
            otherwise
                error('face3DShellDKTTri3N: error face index');
        end
        
        % project to local coordinate
        vector_pressure_local=vector_pressure*lamada';
        
        if EDGE_FLAG
            % edge load
            equivalent_force=zeros(1,DOF_element);
            
            % edge load
            delta_x=x_local(face_node_index(2))-x_local(face_node_index(1));
            delta_y=y_local(face_node_index(2))-y_local(face_node_index(1));
            
            len_edge=norm([delta_x,delta_y]);
            
            b=[y_local(2)-y_local(3),y_local(3)-y_local(1),y_local(1)-y_local(2)];
            c=[x_local(3)-x_local(2),x_local(1)-x_local(3),x_local(2)-x_local(1)];
            
            equivalent_force(6*face_node_index(1)-5:6*face_node_index(1))=vector_pressure_local*...
                [len_edge/2 0 0 0 0 b(bc_index)/12*len_edge;
                0 len_edge/2 0 0 0 c(bc_index)/12*len_edge;
                0 0 len_edge/2 0 0 0]*t_local(face_node_index(1));
            
            equivalent_force(6*face_node_index(2)-5:6*face_node_index(2))=vector_pressure_local*...
                [len_edge/2 0 0 0 0 -b(bc_index)/12*len_edge;
                0 len_edge/2 0 0 0 -c(bc_index)/12*len_edge;
                0 0 len_edge/2 0 0 0]*t_local(face_node_index(2));
            
            % project to global coordinate
            for node_index=1:node_number
                equivalent_force((node_index-1)*DOF_node+1:node_index*DOF_node)=...
                    equivalent_force((node_index-1)*DOF_node+1:node_index*DOF_node)*...
                    [lamada,zeros(3);zeros(3),lamada];
            end
            
            equivalent_force=reshape(equivalent_force,DOF_node,node_number)';
        else
            % surface load
            equivalent_force=zeros(1,DOF_element);
            % stress integral
            for integral_index=1:size(x2_list,1)
                x=x2_list(integral_index,:)*x_local;
                y=x2_list(integral_index,:)*y_local;
                [N,A]=calMatrixN...
                    (element_index,x,y,zeta);
                equivalent_force=equivalent_force+...
                    w2_list(integral_index)*vector_pressure_local*N*A;
            end
            % project to global coordinate
            for node_index=1:node_number
                equivalent_force((node_index-1)*DOF_node+1:node_index*DOF_node)=...
                    equivalent_force((node_index-1)*DOF_node+1:node_index*DOF_node)*...
                    [lamada,zeros(3);zeros(3),lamada];
            end
            equivalent_force=reshape(equivalent_force,DOF_node,node_number)';
        end
    otherwise
        error('face3DShellDKTTri3N: undefined face load');
end
end
function equivalent_force=volume3DShellDKTTri3N(parameter)
% calculate the equivalent nodal force of the distributed pressure
% notice that only normal pressure is available
% use gaussian integal
%
% return global coordinate equivalent_force(node_number*DOF_node matrix)
%
% equivalent_force format:
% {F_x1,F_y1,F_z1,M_x1,M_y1,M_y1;...}
%
% quadratic Gaussian integral point and weight
global g_Element g_Local g_Point_Local

% one Gaussian integral point and weight
x1_list=[0.3333333333333333,0.3333333333333333,0.3333333333333333];
w1_list=1;
% quadratic Gaussian integral point and weight
% x2_list=[
%     0.3333333333333333,0.3333333333333333,0.3333333333333333;
%     0.0000000000000000,0.5000000000000000,0.5000000000000000;
%     0.5000000000000000,0.5000000000000000,0.0000000000000000;
%     0.5000000000000000,0.0000000000000000,0.5000000000000000;];
% w2_list=[0.5;0.1666666666666667;0.1666666666666667;0.1666666666666667];

x2_list=[
    0.6666666666666667,0.1666666666666667,0.1666666666666667;
    0.1666666666666667,0.6666666666666667,0.1666666666666667;
    0.1666666666666667,0.1666666666666667,0.6666666666666667;];
w2_list=[0.3333333333333333;0.3333333333333333;0.3333333333333333];

% quintic Simpson
x5s_list=[-1,              -0.5,            0,               0.5,             1];
w5s_list=[0.166666666666667,0.666666666666667,0.333333333333333,0.666666666666667,0.166666666666667];

node_number=3;
DOF_node=6;
DOF_element=DOF_node*node_number;

element_index=parameter(1);
acceleration=parameter(2:4);

% obtain face normal vector
vector_x=g_Local(element_index*3-2,:);
vector_y=g_Local(element_index*3-1,:);
vector_z=g_Local(element_index*3,:);

lamada=[vector_x;vector_y;vector_z];

x_local=g_Point_Local(node_number*(element_index-1)+1:node_number*element_index,1);
y_local=g_Point_Local(node_number*(element_index-1)+1:node_number*element_index,2);
t_local=g_Point_Local(node_number*(element_index-1)+1:node_number*element_index,3);

% project to local
acceleration_local=acceleration*lamada';

equivalent_force=zeros(1,DOF_element);
for k=1:length(x5s_list)
    for integral_index=1:size(x2_list,1)
        x=x2_list(integral_index,:)*x_local;
        y=x2_list(integral_index,:)*y_local;
        [N,A]=calMatrixN...
            (element_index,x,y,x5s_list(k));
        equivalent_force=equivalent_force+...
            w2_list(integral_index)*w5s_list(k)*acceleration_local*N*A;
    end
end
% project to global coordinate
for node_index=1:node_number
    equivalent_force((node_index-1)*DOF_node+1:node_index*DOF_node)=...
        equivalent_force((node_index-1)*DOF_node+1:node_index*DOF_node)*...
        [lamada,zeros(3);zeros(3),lamada];
end
equivalent_force=reshape(equivalent_force,DOF_node,node_number)';

end
function element_stress=stress3DShellDKTTri3N(element_index,reduced_integrate_flag)
% calculate element stress
% point_index_of_element,{stress_xx,stress_yy,stress_zz(0),stress_xy,stress_yz,stress_zx,stress_von}
%
global g_Element g_Local g_Point_Local g_Delta
node_number=3;
DOF_node=6;
DOF_element=node_number*DOF_node;

node_index=g_Element(element_index,2:4);

% obtain face vector
vector_x=g_Local(element_index*3-2,:);
vector_y=g_Local(element_index*3-1,:);
vector_z=g_Local(element_index*3,:);

% rank_vector*lamada' means global to local
% lamada*colume_vector means global to local
lamada=[vector_x;vector_y;vector_z];

delta=zeros(node_number*DOF_node,1);
delta(1:DOF_node:end)=g_Delta((node_index-1)*DOF_node+1);
delta(2:DOF_node:end)=g_Delta((node_index-1)*DOF_node+2);
delta(3:DOF_node:end)=g_Delta((node_index-1)*DOF_node+3);
delta(4:DOF_node:end)=g_Delta((node_index-1)*DOF_node+4);
delta(5:DOF_node:end)=g_Delta((node_index-1)*DOF_node+5);
delta(6:DOF_node:end)=g_Delta((node_index-1)*DOF_node+6);

x_local=g_Point_Local(node_number*(element_index-1)+1:node_number*element_index,1);
y_local=g_Point_Local(node_number*(element_index-1)+1:node_number*element_index,2);
t_local=g_Point_Local(node_number*(element_index-1)+1:node_number*element_index,3);

% project delta to delta_local
delta_local=lamada*reshape(delta,3,DOF_element/3);
delta_local=delta_local(:);

D=calMatrixD(element_index);

element_stress=zeros(node_number,7);
if reduced_integrate_flag
    B=calMatrixB(element_index,...
        sum(x_local)/node_number,sum(y_local)/node_number,1);
    
    % calculate element stress of up surface
    sigma_local=D*B*delta_local;
    sigma_local=[sigma_local(1:2);0;sigma_local(3);0;0];
    sigma_matrix_local=[
        sigma_local(1) sigma_local(4) sigma_local(5);
        sigma_local(4) sigma_local(2) sigma_local(6);
        sigma_local(5) sigma_local(6) sigma_local(3);];
    
    % principal main stress
    stress_principal=eig(sigma_matrix_local);
    von_sigma=sqrt((...
        (stress_principal(1)-stress_principal(2))^2+...
        (stress_principal(2)-stress_principal(3))^2+...
        (stress_principal(3)-stress_principal(1))^2)/2);
    
    % project sigma_matrix_local to global
    sigma_matrix=lamada'*sigma_matrix_local;
    sigma=[sigma_matrix(1,1),sigma_matrix(2,2),sigma_matrix(3,3),sigma_matrix(1,2),sigma_matrix(1,3),sigma_matrix(2,3)];
    
    element_stress(:,:)=repmat([sigma,von_sigma],node_number,1);
else
    for node_index=1:node_number
        B=calMatrixB(element_index,...
            x_local(node_index),y_local(node_index),1);
        
        % calculate element stress of up surface
        sigma_bend_local=D*B*delta_local;
        sigma_local=[sigma_bend_local(1:2);0;sigma_bend_local(3);0;0];
%         sigma_local=[sigma_bend_local(1:2);0;sigma_bend_local(3);sigma_shear_local];
        sigma_matrix_local=[
            sigma_local(1) sigma_local(4) sigma_local(5);
            sigma_local(4) sigma_local(2) sigma_local(6);
            sigma_local(5) sigma_local(6) sigma_local(3);];
        
        % principal main stress
        stress_principal=eig(sigma_matrix_local);
        von_sigma=sqrt((...
            (stress_principal(1)-stress_principal(2))^2+...
            (stress_principal(2)-stress_principal(3))^2+...
            (stress_principal(3)-stress_principal(1))^2)/2);
        
        % project sigma_matrix_local to sigma_matrix
        sigma_matrix=lamada'*sigma_matrix_local;
        sigma=[sigma_matrix(1,1),sigma_matrix(2,2),sigma_matrix(3,3),sigma_matrix(1,2),sigma_matrix(1,3),sigma_matrix(2,3)];
        
        element_stress(node_index,:)=[sigma,von_sigma];
    end
end
end