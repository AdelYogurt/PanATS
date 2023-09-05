function [max_heat_flux]=solveModelHypersonicHeat()
% plate reference enthalpy method to calculate heat density of surface element
%
% point_list is coordinate of all node
% element_list contain []
% marker_list contain make{marker_name,marker_element_number,marker_element}
% marker_element contain HATSElement
%
% note P_e == P_w because of
% pressure of outer boundary of boundary layer is equal to pressure of surface
%
% copyright Adel 2023.03
%
global user_model

geometry=user_model.geometry;

dimension=geometry.dimension;
point_list=geometry.point_list;
element_list=user_model.element_list;

SYMMETRY=user_model.SYMMETRY;

% heat calculate need inviscid and streamline result
output_inviscid=user_model.output_inviscid;
output_streamline=user_model.output_streamline;
output_boulay=user_model.output_boulay;
output_viscid=user_model.output_viscid;

% calculate inflow vector
free_flow_vector=calFreeFlowDirection(user_model.AOA,user_model.SIDESLIP_ANGLE);
user_model.free_flow_vector=free_flow_vector;

% reference value
T_1=user_model.FREESTREAM_TEMPERATURE;
P_1=user_model.FREESTREAM_PRESSURE;
Ma_1=user_model.MACH_NUMBER;
gamma=user_model.GAMMA_VALUE;
Re=user_model.REYNOLDS_NUMBER;
T_w=user_model.MARKER_ISOTHERMAL;

% load data from inviscid, streamline, boundary layer and viscid result
theta_list=output_inviscid.theta_list;
P_list=output_inviscid.P_list;

boolean_attach_list=output_streamline.boolean_attach_list;
element_attach_list=output_streamline.element_attach_list;

T_2_list=output_boulay.T_2_list;
P_2_list=output_boulay.P_2_list;
rho_2_list=output_boulay.rho_2_list;

rho_e_list=output_boulay.rho_e_list;
V_e_list=output_boulay.V_e_list;
mu_e_list=output_boulay.mu_e_list;
H_e_list=output_boulay.H_e_list;

rho_ref_list=output_boulay.rho_ref_list;
mu_ref_list=output_boulay.mu_ref_list;
H_ref_list=output_boulay.H_ref_list;

H_r_list=output_boulay.H_r_list;

Re_x_list=output_boulay.Re_x_list;
Re_x_ref_list=output_boulay.Re_x_ref_list;

Cf_list=output_viscid.Cf_list;

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
rho_1=P_1/R/T_1;
a_1=sqrt(gamma*R*T_1);
V_1=a_1*Ma_1;
V_1_sq=V_1*V_1;
q_1=rho_1*V_1_sq/2;

% Sutherland method
mu_1=airViscosity(airEnthalpy(T_1));
mu_w=airViscosity(airEnthalpy(T_w));

% solve prepare
H_w=airEnthalpy(T_w); % wall air enthalpy J/kg
H_D=(33867*0.78+15320*0.22)*1e3; % Mean dissociation enthalpy of air J/kg

max_heat_flux=0;

elem_num=length(element_list);
% initialize result sort array
HF_list=zeros(elem_num,1);

% calculate heat flow by plate reference enthalpy method
elem_num=length(element_list);
for elem_idx=1:elem_num
    % if element is cross by attachment line, do not need to calculate
    if boolean_attach_list(elem_idx)
        continue;
    end

    % load data
    Cf=Cf_list(elem_idx);
    H_r=H_r_list(elem_idx);
    rho_e=rho_e_list(elem_idx);
    V_e=V_e_list(elem_idx);

    HF=0.5*Cf*Prpn23*rho_e*V_e*(H_r-H_w);

    HF_list(elem_idx)=HF;
end

% calculate attachment element by Detra-Kemp-Riddell method
for attach_idx=1:length(element_attach_list)
    elem_idx=element_attach_list(attach_idx);
    elem=element_list(elem_idx);

    % stagnation point air parameter parpare
    P_es=P_list(elem_idx); % stagnation point pressure

    T_2s=T_2_list(elem_idx);
    P_2s=P_2_list(elem_idx);
    rho_2s=rho_2_list(elem_idx);

    rho_es=rho_e_list(elem_idx);
    V_es=V_e_list(elem_idx);
    mu_es=mu_e_list(elem_idx);
    H_es=H_e_list(elem_idx);

    P_ws=P_es;
    rho_ws=airDensity(H_w,P_ws);

    % load nearby point coordination to fit ellipsoid
    vertex_arou_idx=[];
    node_num=elem.node_num;
    for node_idx=1:node_num
        [~,Vertex_arou]=getAroundElement(element_list,elem_idx,node_idx);
        vertex_arou_idx=[vertex_arou_idx,Vertex_arou];
    end
    vertex_arou_idx=unique(vertex_arou_idx);
    point_arou_list=point_list(vertex_arou_idx,:);

    % fitting ellipsoid
    %     [~,radius_s_list]=fitEllipsoid(point_arho_list);
    %     radius_s=mean(radius_s_list);radius_s_min=min(radius_s_list);radius_s_max=max(radius_s_list);
    % fitting sphere or cylinder
    [~,radius_shpere,error_shpere]=fitSphere(point_arou_list);
    [direction,radius_cylinder,~,error_cylinder]=fitCylinder(point_arou_list);

    %     HF_s=calFayRiddell(radius_s);
    %     HF_s=calTauber(radius_s);
    %     HF_s=calKempRiddell(radius_s);

    %     if any(sum(abs(point_list(Node_idx,1:3)-[3.175,2.35663,0]),2) < 1e-3)
    %         disp('?');
    %     end

    if error_shpere < error_cylinder
        % use shpere
        HF_s=calDetraKempRiddell(radius_shpere,radius_shpere);
    else
        % use cylinder
        HF_s=calDetraKempRiddell(radius_cylinder,radius_cylinder);
        sin_sweep=abs(dot(direction,free_flow_vector));
        sin_sweep_sq=sin_sweep*sin_sweep;
        HF_s=HF_s*sqrt((1-sin_sweep_sq)^1.5/2);
    end

    %     % stagnation parameter
    %     T_0=T_1*(gamma_sub/2*Ma_1*Ma_1+1);
    %     P_0=P_1*(gamma_sub/2*Ma_1*Ma_1+1)^(gamma/gamma_sub);
    %     rho_0=rho_1*(gamma_sub/2*Ma_1*Ma_1+1)^(1/gamma_sub);
    %     mu_0=airViscosity(airEnthalpy(T_0));
    %     T_es=T_1*(gamma_sub/2*Ma_1*Ma_1+1);
    %     rho_es=rho_1*(gamma_plus*Ma_1*Ma_1)/(gamma_sub*Ma_1*Ma_1+2);
    %     P_es=P_1*((gamma_plus^gamma_plus*Ma_1^(2*gamma))/...
    %         (2^gamma_plus*gamma*Ma_1*Ma_1-gamma_sub))^(1/gamma_sub);
    %     P_w=P_1*(2*gamma*Ma_1*Ma_1-gamma_sub)/gamma_plus;
    %     rho_w=rho_es;
    %     mu_es=airViscosity(airEnthalpy(T_es));

    %     HF_s=calFayRiddell(radius_s_min);

    %     % load arhond element
    %     V_e_aruo_list=[];
    %     center_point_arho_list=[];
    %     Node_idx=element.Node_idx;
    %     point_number=length(Node_idx);
    %     for point_index=1:point_number
    %         vertex_index=Node_idx(point_index);
    %         if point_index == point_number
    %             vertex_ref_index=Node_idx(1);
    %         else
    %             vertex_ref_index=Node_idx(point_index+1);
    %         end
    %
    %         edge_arho=edge_list(vertex_ref_index);
    %         index=edge_arho.getRefIndex(vertex_index);
    %
    %         if ~isempty(index)
    %             element_arho=edge_arho.element_ref_list(index);
    %
    %             V_e_aruo_list=[V_e_list{element_arho.marker_index}(element_arho.element_index);V_e_aruo_list];
    %             center_point_arho_list=[element_arho.center_point;center_point_arho_list];
    %         end
    %     end
    %
    %     % calculate du_e__ds by arhond V_e
    %     dist_list=sqrt(sum((center_point_arho_list-element.center_point).^2,2));
    %     du_e__ds=sum(V_e_aruo_list./dist_list)/point_number;
    %     HF_s=0.763*Pr^-0.6*((rho_ws*mu_w)/(rho_es*mu_es))^0.1*...
    %         sqrt(rho_es*mu_es*du_e__ds)*(1+(Le^a-1)*H_D/H_es)*(H_es-H_w);

    HF_list(elem_idx)=HF_s;

    if max_heat_flux < HF_s
        max_heat_flux=HF_s;
    end

end

output_heat.HF_list=HF_list;

user_model.output_heat=output_heat;

if user_model.INFORMATION
    fprintf('solveModelHypersonicHeat: hypersonic heat solve done!\n');
    fprintf('solveModelHypersonicHeat: result\n');
    fprintf('max heat flux: %14f\n',max_heat_flux)
end

    function HF_s=calFayRiddell(radius_s)
        % Fay-Riddell method
        du_e__ds=sqrt(2*(P_es-P_1)/rho_es)/radius_s;
        HF_s=0.763*Pr^-0.6*((rho_ws*mu_w)/(rho_es*mu_es))^0.1*...
            sqrt(rho_es*mu_es*du_e__ds)*(1+(Le^a-1)*H_D/H_es)*(H_es-H_w);
    end

    function HF_s=calTauber(radius_s)
        % Tauber method
        HF_s=1.83e-4/sqrt(radius_s)*sqrt(rho_1)*V_1^3*(1-H_w/H_es);
    end

    function HF_s=calKempRiddell(radius_s)
        % Kemp-Riddell method
        % HF_s=1.103317e8/sqrt(radius_s)*...
        %     (rho_1/rho_sl)^0.5*(V_1/7925)^3.15*(H_es-H_w)/(H_es-3.0145e5);
        % modify
        HF_s=1.103317e8/sqrt(radius_s)*sqrt(2)*...
            (rho_1/rho_sl)^0.5*(V_1/7925)^3.5*(H_es-H_w)/(H_es-3.0145e5);
    end

    function HF_s=calDetraKempRiddell(radius_s_min,radius_s_max)
        % Detra-Kemp-Riddell method
        HF_s=1.1037e8/sqrt(radius_s_min)*sqrt(1.1+0.9*sqrt(radius_s_min/radius_s_max))*...
            (rho_1/rho_sl)^0.5*(V_1/7925)^3.5*(H_es-H_w)/(H_es-3.0145e5);
    end
end

