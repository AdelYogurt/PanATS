function [Cl,Cd,LDratio,Cx,Cy,Cz,Cmx,Cmy,Cmz]=solveModelHypersonicViscid()
% plate reference enthalpy method to calculate viscid
%
% copyright Adel 2023.03
%
global user_model

geometry_torlance=1e-9;

geometry=user_model.geometry;

dimension=geometry.dimension;
point_list=geometry.point_list;
element_list=user_model.element_list;

SYMMETRY=user_model.SYMMETRY;

% load geometry
center_point_list=geometry.center_point_list;
area_list=geometry.area_list;
normal_vector_list=geometry.normal_vector_list;

% heat calculate need inviscid and streamline result
output_inviscid=user_model.output_inviscid;
output_streamline=user_model.output_streamline;
output_boulay=user_model.output_boulay;

% calculate inflow vector
free_flow_vector=calFreeFlowDirection(user_model.AOA,user_model.SIDESLIP_ANGLE);
user_model.free_flow_vector=free_flow_vector;

% reference value
ref_point=[user_model.REF_ORIGIN_MOMENT_X,user_model.REF_ORIGIN_MOMENT_Y,user_model.REF_ORIGIN_MOMENT_Z];
ref_area=user_model.REF_AREA;
ref_length=user_model.REF_LENGTH;

T_1=user_model.FREESTREAM_TEMPERATURE;
P_1=user_model.FREESTREAM_PRESSURE;
Ma_1=user_model.MACH_NUMBER;
gamma=user_model.GAMMA_VALUE;
Re=user_model.REYNOLDS_NUMBER;

% load data from inviscid and boundary layer result
force_inviscid=output_inviscid.force_inviscid;
moment_inviscid=output_inviscid.moment_inviscid;

surface_flow_list=output_streamline.surface_flow_list;

rho_e_list=output_boulay.rho_e_list;
rho_ref_list=output_boulay.rho_ref_list;

Re_x_list=output_boulay.Re_x_list;
Re_x_ref_list=output_boulay.Re_x_ref_list;

% air parameter
rho_sl=1.225;

R=287.0955;

% free flow parameter
rho_1=P_1/R/T_1;
a_1=sqrt(gamma*R*T_1);
V_1=a_1*Ma_1;
V_1_sq=V_1*V_1;
q_1=rho_1*V_1_sq/2;

% solve prepare
Re_x_tri=10^(5.37+0.2326*Ma_1-0.004015*Ma_1*Ma_1); % transition Reynolds number

elem_num=length(element_list);
% initialize result sort array
Cf_list=zeros(elem_num,1);
S_list=zeros(elem_num,1);
dFs_list=zeros(elem_num,3);
dMs_list=zeros(elem_num,3);

elem_num=length(element_list);
for elem_idx=1:elem_num
    center_point=center_point_list(elem_idx,:);
    area=area_list(elem_idx);
    surface_flow=surface_flow_list(elem_idx,:);

    % load data
    rho_e=rho_e_list(elem_idx);
    rho_ref=rho_ref_list(elem_idx);

    Re_x=Re_x_list(elem_idx);
    Re_x_ref=Re_x_ref_list(elem_idx);

    if Re_x <= Re_x_tri
        % laminar flow
        Cf_ref=0.6640*Re_x_ref^-0.5;
    elseif Re_x < 1e7
        % transition flow
        Cf_ref=0.0296*Re_x_ref^-0.2;
    else
        % turbulent flow
        Cf_ref=0.288*(log(Re_x_ref)/log(10))^2.584;
    end
    Cf=Cf_ref*rho_ref/rho_e;

    % Shear stress
    S=q_1*Cf;

    norm_surface_flow=norm(surface_flow);

    if norm_surface_flow < geometry_torlance
        Cf_list(elem_idx,:)=0;
        S_list(elem_idx,:)=0;
        dFs_list(elem_idx,:)=zeros(1,3);
        dMs_list(elem_idx,:)=zeros(1,3);
    else
        Cf_list(elem_idx,:)=Cf;
        S_list(elem_idx,:)=S;
        dFs_list(elem_idx,:)=surface_flow/norm_surface_flow*area*S;
        dMs_list(elem_idx,:)=cross(center_point-ref_point,dFs_list(elem_idx,:));
    end
end
force_viscid=sum(dFs_list,1);
moment_viscid=sum(dMs_list,1);

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