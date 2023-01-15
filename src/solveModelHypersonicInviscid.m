function [Cl,Cd,LDratio,Cmy]=solveModelHypersonicInviscid...
    (rou_1,V_1,T_1,P_1,Ma_1,gama,AOA,SIDESLIP_ANGLE,...
    ref_point,ref_length,ref_area,Re)
% Newton method to calculate hypersonic aircraft
%
% point_list is coordinate of all node
% element_list contain element(element_type, node_index1, node_index2, ...)
% marker_list contain make{marker_name,marker_element_number,marker_element}
% marker_element contain contain element(element_type, node_index1, node_index2, ...)
%
% copyright Adel 2022.11
%
global g_geometry g_Point g_Element g_Marker...
    ADtree_marker_element HATS_element_list...
    streamline_output inviscid_output

marker_element_number=size(HATS_element_list,1);
dimension=3;

vector_flow=g_geometry.vector_flow;

R=287.0955;
rou_1=P_1/R/T_1;
a_1=sqrt(gama*R*T_1);
V_1=a_1*Ma_1;
q_1=rou_1*V_1*V_1/2;

g_geometry.rou_1=rou_1;
g_geometry.V_1=V_1;
g_geometry.T_1=T_1;
g_geometry.P_1=P_1;
g_geometry.Ma_1=Ma_1;
g_geometry.gama=gama;
g_geometry.Re=Re;

% solve prepare
Ma_1_sq=Ma_1*Ma_1;
gama_plus=gama+1;
gama_sub=gama-1;

% modified Dahlem Buck parameter
a=(6.0-0.3*Ma_1)+sin((log(Ma_1)-0.588)/1.2*pi);
n=-1.15-0.5*sin((log(Ma_1)-0.916)/3.29*pi);

% Dejarnetle parameter
G=2.6054749-0.465998*Ma_1+0.09309305*Ma_1*Ma_1+...
    -0.00817329*Ma_1*Ma_1*Ma_1+0.00026447*Ma_1*Ma_1*Ma_1*Ma_1;
D=1.0081057-0.0132323*Ma_1+0.00164956*Ma_1*Ma_1-0.00006797*Ma_1*Ma_1*Ma_1;

% Prandtl Mayer parameter
delta_max=(sqrt(gama_plus/gama_sub)-1)*pi/2-...
    sqrt(gama_plus/gama_sub)*atan(sqrt(gama_sub/gama_plus*(Ma_1_sq-1)))+...
    atan(sqrt(Ma_1_sq-1));

% modified Newton Cp_max
Cp_max=2/gama/Ma_1_sq*(...
    (gama_plus^2*Ma_1_sq/(4*gama*Ma_1_sq-2*gama_sub))^(gama/gama_sub)*((1-gama+2*gama*Ma_1_sq)/gama_plus)...
    -1);

if Ma_1 < 5
    warning('solveModelHypersonicAerodynamic: Mach number less than 5');
elseif Ma_1 < 8.5
    HIGH_HYPERSONIC_FLAG=0;
else
    HIGH_HYPERSONIC_FLAG=1;
end

delta_list=zeros(marker_element_number,1);
Cp_list=zeros(marker_element_number,1);
P_list=zeros(marker_element_number,1);
dFn_list=zeros(marker_element_number,3);
dMn_list=zeros(marker_element_number,3);
for element_index=1:marker_element_number
    center_point=ADtree_marker_element.center_point_list(element_index,:);
    
    normal_vector=HATS_element_list(element_index).normal_vector;
    area=HATS_element_list(element_index).area;
    
    % theta is angle of negative normal vector and nomlz_vec_air
    cos_theta=-normal_vector*vector_flow;
    if cos_theta < -1
        theta=pi;
    elseif cos_theta > 1
        theta=0;
    else
        theta=acos(cos_theta);
    end
    % delta is attack angle
    delta=abs(pi/2-theta);
    
    % Cp
    if theta < pi/2
        % modified Newton
        Cp=Cp_max*cos_theta^2;
    else
        % ACM experantial
        Cp=max(-(delta*3.8197)*(1/Ma_1_sq),-(1/Ma_1_sq));
    end
    
    % p
    P=q_1*Cp+P_1;
    
    if P < 0
        disp('P < 0');
    end
    
    delta_list(element_index,:)=delta;
    Cp_list(element_index,:)=Cp;
    P_list(element_index,:)=P;
    dFn_list(element_index,:)=-normal_vector*area*P;
    dMn_list(element_index,:)=cross(center_point-ref_point,dFn_list(element_index,:));
end

% calculate lift and drag
data_force=sum(dFn_list,1);
force=sum(data_force,1);
drag=force*vector_flow;
rotation_matrix=[0,0,1;
    0,1,0;
    -1,0,0]';
lift=force*(rotation_matrix*vector_flow);
Cl=lift/ref_area/q_1;
Cd=drag/ref_area/q_1;
LDratio=Cl/Cd;

% calculate moment
moment=sum(dMn_list,1);
moment_y=moment*[0;1;0];
Cmy=moment_y/ref_area/ref_length/q_1;

inviscid_output.delta_list=delta_list;
inviscid_output.Cp_list=Cp_list;
inviscid_output.P_list=P_list;
inviscid_output.dFn_list=dFn_list;
inviscid_output.dMn_list=dMn_list;
inviscid_output.Cl=Cl;
inviscid_output.Cd=Cd;
inviscid_output.Cmz=Cmy;
end