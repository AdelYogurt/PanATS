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

T_inf=user_model.FREESTREAM_TEMPERATURE;
P_inf=user_model.FREESTREAM_PRESSURE;
Ma_inf=user_model.MACH_NUMBER;
gamma=user_model.GAMMA_VALUE;
Re_inf=user_model.REYNOLDS_NUMBER;

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
rho_inf=P_inf/R/T_inf;
a_inf=sqrt(gamma*R*T_inf);
V_inf=a_inf*Ma_inf;
V_inf_sq=V_inf*V_inf;
q_inf=rho_inf*V_inf_sq/2;

% solve prepare
Re_x_tri=10^(5.37+0.2326*Ma_inf-0.004015*Ma_inf*Ma_inf); % transition Reynolds number

% calculate frictional drag coefficient
Cf_ref_list=zeros(size(Re_x_ref_list));

% laminar flow
idx=find(Re_x_list < Re_x_tri);
Cf_ref_list(idx)=0.6640*Re_x_ref_list(idx).^-0.5;

% transition flow
idx=find(Re_x_tri <= Re_x_list & Re_x_list < 1e7);
Cf_ref_list(idx)=0.0592*Re_x_ref_list(idx).^-0.2;

% turbulent flow
idx=find(1e7 <= Re_x_list);
Cf_ref_list(idx)=0.37*(log(Re_x_ref_list(idx))./log(10)).^-2.584;

Cf_list=Cf_ref_list.*rho_ref_list./rho_e_list;

% shear stress (Pa)
S_list=q_inf*Cf_list;

% flow vector
tangent_vector_list=surface_flow_list./vecnorm(surface_flow_list,2,2);
tangent_vector_list(isnan(tangent_vector_list))=0;

% shear stress vector (N)
dFs_list=tangent_vector_list.*area_list.*S_list;
dMs_list=cross(center_point_list-ref_point,dFs_list);

% total force
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
Cl=lift/ref_area/q_inf;
Cd=drag/ref_area/q_inf;
LDratio=Cl/Cd;

% calculate force coefficient
Cx=force(1)/ref_area/q_inf;
Cy=force(2)/ref_area/q_inf;
Cz=force(3)/ref_area/q_inf;

% calculate moment
Cmx=moment(1)/ref_area/ref_length/q_inf;
Cmy=moment(2)/ref_area/ref_length/q_inf;
Cmz=moment(3)/ref_area/ref_length/q_inf;

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