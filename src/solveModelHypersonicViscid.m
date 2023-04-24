function [Cl,Cd,LDratio,Cx,Cy,Cz,Cmx,Cmy,Cmz]=solveModelHypersonicViscid()
% plate reference enthalpy method to calculate viscid
%
% point_list is coordinate of all node
% element_list contain []
% marker_list contain make{marker_name,marker_element_number,marker_element}
% marker_element contain HATSElement
%
% copyright Adel 2023.03
%
global user_model

geometry_torlance=1e-9;

dimension=user_model.geometry.dimension;

point_list=user_model.point_list;
edge_list=user_model.edge_list;
marker_list=user_model.marker_list;

MARKER_MONITORING=user_model.MARKER_MONITORING;
SYMMETRY=user_model.SYMMETRY;

% heat calculate need inviscid and streamline result
output_inviscid=user_model.output_inviscid;
output_streamline=user_model.output_streamline;
output_heat=user_model.output_heat;

% calculate inflow vector
free_flow_vector=[1;0;0];

AOA=user_model.AOA/180*pi;
cos_AOA=cos(AOA);
sin_AOA=sin(AOA);
rotation_AOA=[
    cos_AOA 0 -sin_AOA;
    0 1 0;
    sin_AOA 0 cos_AOA];

AOS=user_model.SIDESLIP_ANGLE/180*pi;
cos_AOS=cos(AOS);
sin_AOS=sin(AOS);
rotation_SIDESLIP_ANGLE=[
    cos_AOS -sin_AOS 0;
    sin_AOS cos_AOS 0;
    0 0 1];

free_flow_vector=rotation_AOA*rotation_SIDESLIP_ANGLE*free_flow_vector;
user_model.free_flow_vector=free_flow_vector;

% reference value
ref_point=[user_model.REF_ORIGIN_MOMENT_X,user_model.REF_ORIGIN_MOMENT_Y,user_model.REF_ORIGIN_MOMENT_Z];
ref_area=user_model.REF_AREA;
ref_length=user_model.REF_LENGTH;

T_1=user_model.FREESTREAM_TEMPERATURE;
P_1=user_model.FREESTREAM_PRESSURE;
Ma_1=user_model.MACH_NUMBER;
gama=user_model.GAMA_VALUE;
Re=user_model.REYNOLDS_NUMBER;
T_w=user_model.MARKER_ISOTHERMAL;

% load data from inviscid and streamline result
Re_x_list=output_heat.Re_x_list;
force_inviscid=output_inviscid.force_inviscid;
moment_inviscid=output_inviscid.moment_inviscid;

% air parameter
rou_sl=1.225;

R=287.0955;
Pr=0.71;
Le=1.4;
a=0.52;
r_l=Pr^(1/3); % temperature recovery coefficient
gama_plus=gama+1;
gama_sub=gama-1;

% free flow parameter
rou_1=P_1/R/T_1;
a_1=sqrt(gama*R*T_1);
V_1=a_1*Ma_1;
V_1_sq=V_1*V_1;
q_1=rou_1*V_1_sq/2;

% solve prepare
Re_x_tri=10^(5.37+0.2326*Ma_1-0.004015*Ma_1*Ma_1); % transition Reynolds number

% initialize result sort array
Cf_list=cell(length(marker_list),1);
S_list=cell(length(marker_list),1);
dFs_list=cell(length(marker_list),1);
dMs_list=cell(length(marker_list),1);
force_viscid=zeros(1,3);
moment_viscid=zeros(1,3);

for monitor_index=1:length(MARKER_MONITORING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITORING(monitor_index),marker_list);

    Cf_list_marker=zeros(marker_list(marker_index).element_number,1);
    S_list_marker=zeros(marker_list(marker_index).element_number,1);
    dFs_list_marker=zeros(marker_list(marker_index).element_number,3);
    dMs_list_marker=zeros(marker_list(marker_index).element_number,3);

    for element_index=1:marker_list(marker_index).element_number
        element=marker_element(element_index);

        area=element.area;
        surface_flow=element.surface_flow;

        Re_x=Re_x_list{marker_index}(element_index);

        if Re_x == 0
            Re_x=10;
        end

        if Re_x <= Re_x_tri
            % laminar flow
            Cf=0.6640*Re_x^-0.5;
%         elseif Re_x < 1e7
%             % transition flow
%             Cf=0.0296*Re_x^-0.2;
        else
            % turbulent flow
            Cf=0.288*(log(Re_x)/log(10))^2.584;
        end

        % Shear stress
        S=q_1*Cf;

        if norm(surface_flow) < geometry_torlance
            Cf_list_marker(element_index,:)=0;
            S_list_marker(element_index,:)=0;
            dFs_list_marker(element_index,:)=zeros(1,3);
            dMs_list_marker(element_index,:)=zeros(1,3);
        else
            Cf_list_marker(element_index,:)=Cf;
            S_list_marker(element_index,:)=S;
            dFs_list_marker(element_index,:)=surface_flow/norm(surface_flow)*area*S;
            dMs_list_marker(element_index,:)=cross(element.center_point-ref_point,dFs_list_marker(element_index,:));
        end

        force_viscid=force_viscid+dFs_list_marker(element_index,:);
        moment_viscid=moment_viscid+dMs_list_marker(element_index,:);
    end

    Cf_list{marker_index}=Cf_list_marker;
    S_list{marker_index}=S_list_marker;
    dFs_list{marker_index}=dFs_list_marker;
    dMs_list{marker_index}=dMs_list_marker;
end

force=force_inviscid+force_viscid;
moment=moment_inviscid+moment_viscid;

% calculate lift and drag coefficient
drag=force*free_flow_vector;
rotation_matrix=[0,0,1;
    0,1,0;
    -1,0,0]';
lift=force*(rotation_matrix*free_flow_vector);
Cl=lift/ref_area/q_1;
Cd=drag/ref_area/q_1;
LDratio=Cl/Cd;

% calculate force coefficient
Cx=force(1)/ref_area/q_1;
Cy=force(2)/ref_area/q_1;
Cz=force(3)/ref_area/q_1;

% calculate moment
Cmx=moment(1)/ref_area/ref_length/q_1;
Cmy=moment(2)/ref_area/ref_length/q_1;
Cmz=moment(3)/ref_area/ref_length/q_1;

% process SYMMETRY
if ~isempty(SYMMETRY)
    switch user_model.SYMMETRY
        case 'XOY'
            Cz=0;
            Cmx=0;
            Cmy=0;
        case 'YOZ'
            Cx=0;
            Cmy=0;
            Cmz=0;
        case 'ZOX'
            Cy=0;
            Cmz=0;
            Cmx=0;
        otherwise
            error('solveModelHypersonicInviscid: nuknown SYMMETRY type');
    end
end

output_viscid.Cf_list=Cf_list;
output_viscid.dFs_list=dFs_list;
output_viscid.dMs_list=dMs_list;
output_viscid.force_viscid=force_viscid;
output_viscid.moment_viscid=moment_viscid;
user_model.output_viscid=output_viscid;

if user_model.INFORMATION
    fprintf('solveModelHypersonicViscid: hypersonic viscid solve done!\n');
    fprintf('solveModelHypersonicViscid: result\n');
    fprintf('Cl:  %14f, Cd:  %14f, L/D: %14f\nCx:  %14f, Cy:  %14f, Cz:  %14f\nCmx: %14f, Cmy: %14f, Cmz: %14f\n',...
        Cl,Cd,LDratio,Cx,Cy,Cz,Cmx,Cmy,Cmz)
end

end