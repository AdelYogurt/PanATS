function [Cl,Cd,LDratio,Cmz]=solveModelHypersonicViscid...
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
    streamline_output inviscid_output heat_output viscid_output

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

Re_x_list=heat_output.Re_x_list;
inflow_direction_list=streamline_output.inflow_direction_list;

dFs_list=zeros(marker_element_number,3);
dMs_list=zeros(marker_element_number,3);
for element_index=1:marker_element_number
    center_point=ADtree_marker_element.center_point_list(element_index,:);
    
    area=HATS_element_list(element_index).area;
    
    Re_x=Re_x_list(element_index);
    inflow_direction=inflow_direction_list(element_index,:);
    
    if Re_x == 0
        Re_x=10;
    end
    
    if Re_x < 2540
        % laminar flow
        Cf=0.6640/Re_x^0.5;
    else
        % turbulent flow
        Cf=0.088*(log(Re_x)-2.3686)/(log(Re_x-1.5)^3);
    end
    
%     % modified Schlichting correlation
%     Cf=0.42/((log(Re_x)/log(10))^2.55*(1+0.25*Ma_1_sq)^0.31);
    
    % Shear
    S=q_1*Cf;
    
    dFs_list(element_index,:)=inflow_direction*area*S;
    dMs_list(element_index,:)=cross(center_point-ref_point,dFs_list(element_index,:));
end

% calculate force
data_force_friction=sum(dFs_list,1);
force=sum(data_force_friction,1);
drag=force*vector_flow;
rotation_matrix=[0,0,1;
    0,1,0;
    -1,0,0]';
lift=force*(rotation_matrix*vector_flow);
Cl=lift/ref_area/q_1+inviscid_output.Cl;
Cd=drag/ref_area/q_1+inviscid_output.Cd;
LDratio=Cl/Cd;

% calculate moment
moment=sum(dMs_list,1);
moment_z=moment*[0;0;1];
Cmz=moment_z/ref_area/ref_length/q_1+inviscid_output.Cmz;

viscid_output.dFs_list=dFs_list;
viscid_output.dMs_list=dMs_list;
viscid_output.Cl=Cl;
viscid_output.Cd=Cd;
viscid_output.Cmz=Cmz;
end