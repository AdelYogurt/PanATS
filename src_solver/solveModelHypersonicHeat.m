function [model,CT_out]=solveModelHypersonicHeat(model)
% plate reference enthalpy method to calculate heat density of surface element
%
% notice:
% P_e == P_w because of pressure of outer boundary of boundary layer is equal to pressure of surface
%
% reference:
% [1] 张志成. 高超声速气动热和热防护[M]. 北京：国防工业出版社, 2003.
% [2] Fay J A, Riddell F R. Theory of Stagnation Point Heat Transfer in
% Dissociated Air[J]. Journal of the Aerospace Sciences, 1958, 25(2):
% 73-85.
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
T_w=config.MARKER_ISOTHERMAL;

% load geometry
dim=geometry.dimension;
pnt_list=geometry.point_list;
elem_list=geometry.element_list;
conn_mat=geometry.connectivity_matrix;
cntr_pnt_list=geometry.center_point_list;

% load numerics
theta_list=numerics.theta_list;
P_list=numerics.P_list;
elem_flow_list=numerics.element_flow_list;
bool_att_list=numerics.boolean_attachment_list;
elem_att_list=numerics.element_attachment_list;
T_p_list=numerics.T_p_list;
P_p_list=numerics.P_p_list;
rho_p_list=numerics.rho_p_list;
rho_e_list=numerics.rho_e_list;
V_e_list=numerics.V_e_list;
mu_e_list=numerics.mu_e_list;
H_e_list=numerics.H_e_list;
rho_ref_list=numerics.rho_ref_list;
mu_ref_list=numerics.mu_ref_list;
H_ref_list=numerics.H_ref_list;
H_r_list=numerics.H_r_list;
Re_x_list=numerics.Re_x_list;
Re_x_ref_list=numerics.Re_x_ref_list;
Cf_list=numerics.Cf_list;

% initialize numerics
numerics.HF_list=[];

% initialize output
CT_out.HFmax=[];

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

% air parameter
rho_sl=1.225;
R=287.0955;
Pr=0.71;
Le=1.4;
a=0.52;
r_l=Pr^(1/3); % temperature recovery coefficient
Prpn23=Pr^(-2/3);
gamma_plus=gamma+1;
gamma_sub=gamma-1;

% free flow parameter
rho_inf=P_inf/R/T_inf;
a_inf=sqrt(gamma*R*T_inf);
V_inf=a_inf*Ma_inf;
V_inf_sq=V_inf*V_inf;
q_inf=rho_inf*V_inf_sq/2;

% Sutherland method
mu_inf=airViscosity(airEnthalpy(T_inf));
mu_w=airViscosity(airEnthalpy(T_w));

% solve prepare
H_w=airEnthalpy(T_w); % wall air enthalpy J/kg
% H_D=(33867*0.78+15320*0.22)*1e3; % Mean dissociation enthalpy of air J/kg
H_D=airEnthalpy(T_inf); % base on calculate result, H_D should be enthalpy of air

elem_num=length(elem_list);
% initialize result sort array
HF_list=zeros(elem_num,1);

% calculate heat flow by plate reference enthalpy method
for elem_idx=1:elem_num
    % if element is cross by attachment line, do not need to calculate
%     if boolean_attach_list(elem_idx)
%         continue;
%     end

    % load data
    Cf=Cf_list(elem_idx);
    H_r=H_r_list(elem_idx);
    rho_e=rho_e_list(elem_idx);
    V_e=V_e_list(elem_idx);

    HF=0.5*Cf*Prpn23*rho_e*V_e*(H_r-H_w);

    HF_list(elem_idx)=HF;
end

% calculate attachment element by Detra-Kemp-Riddell method
for attach_idx=1:length(elem_att_list)
    elem_idx=elem_att_list(attach_idx);

    if isnan(elem_list(elem_idx,4)),node_idx_list=elem_list(elem_idx,1:3);node_num=3;
    else,node_idx_list=elem_list(elem_idx,1:4);node_num=4;end

    % load around element properties
    elem_arou_idx_idx=zeros(1,node_num);
    for node_idx=1:node_num
        pnt_idx=node_idx_list(node_idx);

        if node_idx == node_num
            pnt_jdx=node_idx_list(1);
        else
            pnt_jdx=node_idx_list(node_idx+1);
        end

        elem_arou_idx_idx(node_idx)=full(conn_mat(pnt_jdx,pnt_idx));
    end
    elem_arou_idx_idx(elem_arou_idx_idx==0)=[];

    % stagnation point air parameter parpare
    P_e=P_list(elem_idx); % stagnation point pressure

    T_2=T_p_list(elem_idx);
    P_2=P_p_list(elem_idx);
    rho_2=rho_p_list(elem_idx);

    rho_e=rho_e_list(elem_idx);
    V_e=V_e_list(elem_idx);
    mu_e=mu_e_list(elem_idx);
    H_e=H_e_list(elem_idx);

    P_w=P_e;
    rho_w=airDensity(H_w,P_w);

    % % load nearby point coordination to fit ellipsoid
    % vertex_arou_idx=[];
    % node_num=elem.node_num;
    % for node_idx=1:node_num
    %     [~,Vertex_arou]=getAroundElement(element_list,elem_idx,node_idx);
    %     vertex_arou_idx=[vertex_arou_idx,Vertex_arou];
    % end
    % vertex_arou_idx=unique(vertex_arou_idx);
    % point_arou_list=point_list(vertex_arou_idx,:);
    % 
    % % fitting sphere or cylinder
    % [~,radius_shpere,error_shpere]=fitSphere(point_arou_list);
    % [direction,radius_cylinder,~,error_cylinder]=fitCylinder(point_arou_list);
    % 
    % if error_shpere < error_cylinder
    %     % use shpere
    %     HF_s=calDetraKempRiddell(radius_shpere,radius_shpere);
    % else
    %     % use cylinder
    %     HF_s=calDetraKempRiddell(radius_cylinder,radius_cylinder);
    %     sin_sweep=abs(dot(direction,free_flow_vector));
    %     sin_sweep_sq=sin_sweep*sin_sweep;
    %     HF_s=HF_s*sqrt((1-sin_sweep_sq)^1.5/2);
    % end

    % point coordinate and velocity vector
    arou_pnt=cntr_pnt_list(elem_arou_idx_idx,:);
    V_e_arou=elem_flow_list(elem_arou_idx_idx,:)./vecnorm(elem_flow_list(elem_arou_idx_idx,:),2,2);
    V_e_arou=V_e_arou.*V_e_list(elem_arou_idx_idx,:);

    % point coordinate and velocity vector
    base_pnt=cntr_pnt_list(elem_idx,:);
    V_e=elem_flow_list(elem_idx,:)./vecnorm(elem_flow_list(elem_idx,:),2,2);
    V_e=V_e*V_e_list(elem_idx,:);
    
    % judge if element lay in symmetry plane
    % if yes, symmetric velocity consider
    if identify_dim ~= 0
        arou_pnt_sym=base_pnt;
        arou_pnt_sym(:,identify_dim)=-arou_pnt_sym(:,identify_dim);
        arou_pnt=[arou_pnt;arou_pnt_sym];

        V_e_arou_sym=V_e;
        V_e_arou_sym(:,identify_dim)=-V_e_arou_sym(:,identify_dim);
        V_e_arou=[V_e_arou;V_e_arou_sym];
    end

    % calculate du_e__ds by arhond V_e
    dist_list=sqrt(sum((arou_pnt-base_pnt).^2,2));
    % du_e__ds=mean(vecnorm(V_e_arou-V_e,2,2)./dist_list);
    du_e__ds=mean(abs(vecnorm(V_e_arou,2,2)-norm(V_e))./dist_list);

    % [2] Page: 85, for modified Newtonian flow
    % du_e__ds=sqrt(2*(P_e-P_inf)/rho_e)/radius_s;

    % Fay-Riddell method
    % [1] Pages: 66-67
    g=H_w/H_e;
    l=0.216./sqrt(g)-0.01657./g;
    rho_e_mu_e=l*(rho_w*mu_w);

    % HF=0.763*Pr^-0.6*((rho_w*mu_w)/(rho_e*mu_e))^0.1*sqrt(rho_e*mu_e*du_e__ds)*(1+(Le^a-1)*H_D/H_e)*(H_e-H_w);
    HF=0.763*Pr^-0.6*((rho_w*mu_w)/rho_e_mu_e)^0.1*sqrt(rho_e_mu_e*du_e__ds)*(1+(Le^a-1)*H_D/H_e)*(H_e-H_w);

    HF_list(elem_idx)=max(HF_list(elem_idx),HF);
end

HFmax=max(HF_list);

%% sort data

if config.INFORMATION
    fprintf('solveModelHypersonicHeat: hypersonic heat solve done!\n');
    fprintf('solveModelHypersonicHeat: result\n');
    fprintf('max heat flux: %14f\n',HFmax)
end

numerics.HF_list=HF_list;

CT_out.HFmax=HFmax;

model.numerics=numerics;

%% aerothermal engineering estimation function

    function HF_s=calTauber(radius_s)
        % Tauber method
        HF_s=1.83e-4/sqrt(radius_s)*sqrt(rho_inf)*V_inf^3*(1-H_w/H_e);
    end

    function HF_s=calKempRiddell(radius_s)
        % Kemp-Riddell method
        % HF_s=1.103317e8/sqrt(radius_s)*...
        %     (rho_1/rho_sl)^0.5*(V_1/7925)^3.15*(H_es-H_w)/(H_es-3.0145e5);
        % modify
        HF_s=1.103317e8/sqrt(radius_s)*sqrt(2)*...
            (rho_inf/rho_sl)^0.5*(V_inf/7925)^3.5*(H_e-H_w)/(H_e-3.0145e5);
    end

    function HF_s=calDetraKempRiddell(radius_s_min,radius_s_max)
        % Detra-Kemp-Riddell method
        HF_s=1.1037e8/sqrt(radius_s_min)*sqrt(1.1+0.9*sqrt(radius_s_min/radius_s_max))*...
            (rho_inf/rho_sl)^0.5*(V_inf/7925)^3.5*(H_e-H_w)/(H_e-3.0145e5);
    end
end

% % stagnation parameter
% T_0=T_1*(gamma_sub/2*Ma_1*Ma_1+1);
% P_0=P_1*(gamma_sub/2*Ma_1*Ma_1+1)^(gamma/gamma_sub);
% rho_0=rho_1*(gamma_sub/2*Ma_1*Ma_1+1)^(1/gamma_sub);
% mu_0=airViscosity(airEnthalpy(T_0));
% T_es=T_1*(gamma_sub/2*Ma_1*Ma_1+1);
% rho_es=rho_1*(gamma_plus*Ma_1*Ma_1)/(gamma_sub*Ma_1*Ma_1+2);
% P_es=P_1*((gamma_plus^gamma_plus*Ma_1^(2*gamma))/...
%     (2^gamma_plus*gamma*Ma_1*Ma_1-gamma_sub))^(1/gamma_sub);
% P_w=P_1*(2*gamma*Ma_1*Ma_1-gamma_sub)/gamma_plus;
% rho_w=rho_es;
% mu_es=airViscosity(airEnthalpy(T_es));
    