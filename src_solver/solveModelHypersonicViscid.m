function [CL,CD,CSF,CFx,CFy,CFz,CMx,CMy,CMz,CEff]=solveModelHypersonicViscid()
% plate reference enthalpy method to calculate viscid
%
% copyright Adel 2023.03
%
global user_model

config=user_model.config;
geometry=user_model.geometry;
element_list=user_model.element_list;
SYMMETRY=config.SYMMETRY;

% load geometry
dimension=geometry.dimension;
point_list=geometry.point_list;
center_point_list=geometry.center_point_list;
area_list=geometry.area_list;
normal_vector_list=geometry.normal_vector_list;

% heat calculate need inviscid and streamline result
output_inviscid=user_model.output_inviscid;
output_streamline=user_model.output_streamline;
output_boulay=user_model.output_boulay;

% calculate inflow vector
coord_vec=coordVecToOriSU2(config.AOA,config.SIDESLIP_ANGLE,dimension);
free_flow_vector=coord_vec(:,1);

% reference value
ref_point=[config.REF_ORIGIN_MOMENT_X,config.REF_ORIGIN_MOMENT_Y,config.REF_ORIGIN_MOMENT_Z];
ref_area=config.REF_AREA;
ref_length=config.REF_LENGTH;

T_inf=config.FREESTREAM_TEMPERATURE;
P_inf=config.FREESTREAM_PRESSURE;
Ma_inf=config.MACH_NUMBER;
gamma=config.GAMMA_VALUE;
% Re_inf=config.REYNOLDS_NUMBER;

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
srf_flow_norm=vecnorm(surface_flow_list,2,2);
srf_flow_norm(srf_flow_norm == 0)=1;
tangent_vector_list=surface_flow_list./srf_flow_norm;

% shear stress vector (N)
dFs_list=tangent_vector_list.*area_list.*S_list;
dMs_list=cross(center_point_list-ref_point,dFs_list);

% total force
force_viscid=sum(dFs_list,1);
moment_viscid=sum(dMs_list,1);

force=force_inviscid+force_viscid;
moment=moment_inviscid+moment_viscid;

% calculate lift and drag coefficient
drag=force*coord_vec(:,1);
slip=force*coord_vec(:,2);
lift=force*coord_vec(:,3);

% calculate velocity coefficient
CL=lift/ref_area/q_inf;
CD=drag/ref_area/q_inf;
CSF=slip/ref_area/q_inf;

% calculate force coefficient
CFx=force(1)/ref_area/q_inf;
CFy=force(2)/ref_area/q_inf;
CFz=force(3)/ref_area/q_inf;

% calculate moment coefficient
CMx=moment(1)/ref_area/ref_length/q_inf;
CMy=moment(2)/ref_area/ref_length/q_inf;
CMz=moment(3)/ref_area/ref_length/q_inf;

% calculate efficient coefficient
CEff=CL/(CD+eps);

% process SYMMETRY
if ~isempty(SYMMETRY)
    switch config.SYMMETRY
        case 'XOY'
            CFz=0;
            CMx=0;
            CMy=0;
        case 'YOZ'
            CFx=0;
            CMy=0;
            CMz=0;
        case 'ZOX'
            CFy=0;
            CMz=0;
            CMx=0;
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

if config.INFORMATION
    fprintf('solveModelHypersonicViscid: hypersonic viscid solve done!\n');
    fprintf('solveModelHypersonicViscid: result\n');
    fprintf('CL:  %14f, CD:  %14f, CSF: %14f\nCFx:  %14f, CFy:  %14f, CFz:  %14f\nCMx: %14f, CMy: %14f, CMz: %14f\nCEff: %14f\n',...
        [CL,CD,CSF,CFx,CFy,CFz,CMx,CMy,CMz,CEff])
end

end