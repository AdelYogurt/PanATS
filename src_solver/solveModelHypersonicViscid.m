function [model,CA_out]=solveModelHypersonicViscid(model)
% plate reference enthalpy method to calculate viscid
%
% copyright Adel 2023.03
%
config=model.config; % load config
geometry=model.geometry; % load geometry
numerics=model.numerics; % load numerics
INFORMATION=config.INFORMATION; % whether print information

% load config
T_inf=config.FREESTREAM_TEMPERATURE;
P_inf=config.FREESTREAM_PRESSURE;
Ma_inf=config.MACH_NUMBER;
gamma=config.GAMMA_VALUE;
SYMMETRY=config.SYMMETRY;
ref_point=[config.REF_ORIGIN_MOMENT_X,config.REF_ORIGIN_MOMENT_Y,config.REF_ORIGIN_MOMENT_Z];
ref_area=config.REF_AREA;
ref_length=config.REF_LENGTH;

% load geometry
dim=geometry.dimension;
pnt_list=geometry.point_list;
elem_list=geometry.element_list;
cntr_pnt_list=geometry.center_point_list;
E_nmvctr_list=geometry.EN_vector_list;
P_nmvctr_list=geometry.PN_vector_list;
area_list=geometry.area_list;

% load numerics
force_inviscid=numerics.force_inviscid;
moment_inviscid=numerics.moment_inviscid;
elem_flow_list=numerics.element_flow_list;
rho_e_list=numerics.rho_e_list;
rho_ref_list=numerics.rho_ref_list;
Re_x_list=numerics.Re_x_list;
Re_x_ref_list=numerics.Re_x_ref_list;

% initialize numerics
numerics.Cf_list=[];
numerics.dFs_list=[];
numerics.dMs_list=[];
numerics.force_viscid=[];
numerics.moment_viscid=[];

% initialize output
CA_out.CD=[];
CA_out.CL=[];
CA_out.CSF=[];
CA_out.CFx=[];
CA_out.CFy=[];
CA_out.CFz=[];
CA_out.CMx=[];
CA_out.CMy=[];
CA_out.CMz=[];
CA_out.CEff=[];

switch SYMMETRY
    case 'XOY'
        identify_dim=3;
    case 'YOZ'
        identify_dim=1;
    case 'ZOX'
        identify_dim=2;
    otherwise
        identify_dim=0;
end

%% pre process

% calculate element free flow vector
coord_vec=coordVecToOri(config.AOA,config.SIDESLIP_ANGLE,dim);

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

%% calculate surface frictional coefficient base on engineering estimation function 

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
srf_flow_norm=vecnorm(elem_flow_list,2,2);
srf_flow_norm(srf_flow_norm == 0)=1;
tang_vctr_list=elem_flow_list./srf_flow_norm;

% shear stress vector (N)
dFs_list=tang_vctr_list.*area_list.*S_list;
dMs_list=cross(cntr_pnt_list-ref_point,dFs_list);

%% post prcocess

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
if identify_dim
    switch identify_dim
        case 3
            CFz=0;
            CMx=0;
            CMy=0;
        case 1
            CFx=0;
            CMy=0;
            CMz=0;
        case 2
            CFy=0;
            CMz=0;
            CMx=0;
    end
end

%% sort data

if config.INFORMATION
    fprintf('solveModelHypersonicViscid: hypersonic viscid solve done!\n');
    fprintf('solveModelHypersonicViscid: result\n');
    fprintf('CD:  %14f, CL:  %14f, CSF: %14f\nCFx:  %14f, CFy:  %14f, CFz:  %14f\nCMx: %14f, CMy: %14f, CMz: %14f\nCEff: %14f\n',...
        [CD,CL,CSF,CFx,CFy,CFz,CMx,CMy,CMz,CEff])
end

numerics.Cf_list=Cf_list;
numerics.dFs_list=dFs_list;
numerics.dMs_list=dMs_list;
numerics.force_viscid=force_viscid;
numerics.moment_viscid=moment_viscid;

CA_out.CD=CD;
CA_out.CL=CL;
CA_out.CSF=CSF;
CA_out.CFx=CFx;
CA_out.CFy=CFy;
CA_out.CFz=CFz;
CA_out.CMx=CMx;
CA_out.CMy=CMy;
CA_out.CMz=CMz;
CA_out.CEff=CEff;

model.numerics=numerics;
end