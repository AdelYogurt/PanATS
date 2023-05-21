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

dimension=user_model.geometry.dimension;

point_list=user_model.point_list;
edge_list=user_model.edge_list;
marker_list=user_model.marker_list;

MARKER_MONITORING=user_model.MARKER_MONITORING;
SYMMETRY=user_model.SYMMETRY;

% heat calculate need inviscid and streamline result
output_inviscid=user_model.output_inviscid;
output_streamline=user_model.output_streamline;
output_boulay=user_model.output_boulay;
output_viscid=user_model.output_viscid;

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
gamma=user_model.GAMMA_VALUE;
Re=user_model.REYNOLDS_NUMBER;
T_w=user_model.MARKER_ISOTHERMAL;

% load data from inviscid, streamline, boundary layer and viscid result
theta_list=output_inviscid.theta_list;
P_list=output_inviscid.P_list;

attachment_list=output_streamline.attachment_list;

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

% initialize result sort array
Q_list=cell(length(marker_list),1);

% calculate heat flow by plate reference enthalpy method
for monitor_index=1:length(MARKER_MONITORING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITORING(monitor_index),marker_list);
    
    Q_list_marker=zeros(marker_list(marker_index).element_number,1);

    for element_index=1:marker_list(marker_index).element_number
        % if element is cross by attachment line, do not need to calculate
        if marker_element(element_index).attachment
            continue;
        end

        % load data
        Cf=Cf_list{marker_index}(element_index);
        H_r=H_r_list{marker_index}(element_index);
        rho_e=rho_e_list{marker_index}(element_index);
        V_e=V_e_list{marker_index}(element_index);

        Q=0.5*Cf*Prpn23*rho_e*V_e*(H_r-H_w);

        Q_list_marker(element_index)=Q;
    end

    Q_list{marker_index}=Q_list_marker;
end

% calculate attachment element by Detra-Kemp-Riddell method
for attachment_index=1:length(attachment_list)
    marker_index=attachment_list(attachment_index,1);
    element_index=attachment_list(attachment_index,2);
    element=marker_list(marker_index).element_list(element_index);

    % stagnation point air parameter parpare
    P_es=P_list{marker_index}(element_index); % stagnation point pressure

    T_2s=T_2_list{marker_index}(element_index);
    P_2s=P_2_list{marker_index}(element_index);
    rho_2s=rho_2_list{marker_index}(element_index);

    rho_es=rho_e_list{marker_index}(element_index);
    V_es=V_e_list{marker_index}(element_index);
    mu_es=mu_e_list{marker_index}(element_index);
    H_es=H_e_list{marker_index}(element_index);

    P_ws=P_es;
    rho_ws=airDensity(H_w,P_ws);

    % load nearby point coordination to fit ellipsoid
    point_arho_index_list=[];
    point_index_list=element.point_index_list;
    for point_index=1:length(point_index_list)
        vertex_index=point_index_list(point_index);
        edge=edge_list(vertex_index);
        point_arho_index_list=[point_arho_index_list;edge.vertex_ref_list(1:edge.edge_number)];
    end
    point_arho_index_list=unique(point_arho_index_list);
    point_arho_list=point_list(point_arho_index_list,1:3);

    % fitting ellipsoid
    %     [~,radius_s_list]=fitEllipsoid(point_arho_list);
    %     radius_s=mean(radius_s_list);radius_s_min=min(radius_s_list);radius_s_max=max(radius_s_list);
    % fitting sphere or cylinder
    [~,radius_shpere,error_shpere]=fitSphere(point_arho_list);
    [direction,radius_cylinder,~,error_cylinder]=fitCylinder(point_arho_list);

    %     Q_s=calFayRiddell(radius_s);
    %     Q_s=calTauber(radius_s);
    %     Q_s=calKempRiddell(radius_s);
    
%     if any(sum(abs(point_list(point_index_list,1:3)-[3.175,2.35663,0]),2) < 1e-3)
%         disp('?');
%     end
    
    if error_shpere < error_cylinder
        % use shpere
        Q_s=calDetraKempRiddell(radius_shpere,radius_shpere);
    else
        % use cylinder
        Q_s=calDetraKempRiddell(radius_cylinder,radius_cylinder);
        sin_sweep=abs(dot(direction,free_flow_vector));
        sin_sweep_sq=sin_sweep*sin_sweep;
        Q_s=Q_s*sqrt((1-sin_sweep_sq)^1.5/2);
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

%     Q_s=calFayRiddell(radius_s_min);

%     % load arhond element
%     V_e_aruo_list=[];
%     center_point_arho_list=[];
%     point_index_list=element.point_index_list;
%     point_number=length(point_index_list);
%     for point_index=1:point_number
%         vertex_index=point_index_list(point_index);
%         if point_index == point_number
%             vertex_ref_index=point_index_list(1);
%         else
%             vertex_ref_index=point_index_list(point_index+1);
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
%     Q_s=0.763*Pr^-0.6*((rho_ws*mu_w)/(rho_es*mu_es))^0.1*...
%         sqrt(rho_es*mu_es*du_e__ds)*(1+(Le^a-1)*H_D/H_es)*(H_es-H_w);

    Q_list{marker_index}(element_index)=Q_s;

    if max_heat_flux < Q_s
        max_heat_flux=Q_s;
    end

end

output_heat.Q_list=Q_list;

user_model.output_heat=output_heat;

if user_model.INFORMATION
    fprintf('solveModelHypersonicHeat: hypersonic heat solve done!\n');
    fprintf('solveModelHypersonicHeat: result\n');
    fprintf('max heat flux: %14f\n',max_heat_flux)
end

    function Q_s=calFayRiddell(radius_s)
        % Fay-Riddell method
        du_e__ds=sqrt(2*(P_es-P_1)/rho_es)/radius_s;
        Q_s=0.763*Pr^-0.6*((rho_ws*mu_w)/(rho_es*mu_es))^0.1*...
            sqrt(rho_es*mu_es*du_e__ds)*(1+(Le^a-1)*H_D/H_es)*(H_es-H_w);
    end

    function Q_s=calTauber(radius_s)
        % Tauber method
        Q_s=1.83e-4/sqrt(radius_s)*sqrt(rho_1)*V_1^3*(1-H_w/H_es);
    end

    function Q_s=calKempRiddell(radius_s)
        % Kemp-Riddell method
        % Q_s=1.103317e8/sqrt(radius_s)*...
        %     (rho_1/rho_sl)^0.5*(V_1/7925)^3.15*(H_es-H_w)/(H_es-3.0145e5);
        % modify
        Q_s=1.103317e8/sqrt(radius_s)*sqrt(2)*...
            (rho_1/rho_sl)^0.5*(V_1/7925)^3.5*(H_es-H_w)/(H_es-3.0145e5);
    end

    function Q_s=calDetraKempRiddell(radius_s_min,radius_s_max)
        % Detra-Kemp-Riddell method
        Q_s=1.1037e8/sqrt(radius_s_min)*sqrt(1.1+0.9*sqrt(radius_s_min/radius_s_max))*...
            (rho_1/rho_sl)^0.5*(V_1/7925)^3.5*(H_es-H_w)/(H_es-3.0145e5);
    end
end

function [center_point,radius]=fitEllipsoid(point_list)

x=point_list(:,1);
y=point_list(:,2);
z=point_list(:,3);

D=[x.*x,y.*y,z.*z,x.*y,x.*z,y.*z,x,y,z];
P=ones(length(x),1);

C=D\P;

M=[C(1) C(4)/2 C(5)/2;
    C(4)/2 C(2) C(6)/2;
    C(5)/2 C(6)/2 C(3)];

center_point=-0.5*M\[C(7);C(8);C(9)];

S=center_point'*M*center_point+1;
lamada=eig(M);

radius=[sqrt(S/lamada(1));sqrt(S/lamada(2));sqrt(S/lamada(3))];

end

function [center_point,radius,error]=fitSphere(point_list)
% least squares to solve the sphere
%
point_number=size(point_list,1);

matrix=[-2*point_list,ones(point_number,1)];
U=-sum(point_list.^2,2);
coeff=(matrix\U); % xc, yc, zc, D
center_point=coeff(1:3);
% D=xc*xc+yc*yc+zc*zc-r*r
radius=sqrt(sum(center_point.^2)-coeff(4));

% calculate error
error_list=sum((point_list-center_point').^2,2)-radius*radius;
error=sum(error_list.^2,1)/point_number;

end

function [direction,radius,base_point,error]=fitCylinder(point_list)
% least squares to solve the cylinder
%
point_number=size(point_list,1);

normal=pca(point_list);
direction=normal(:,1)';
C=sum(direction.^2);

X=point_list(:,1);
Y=point_list(:,2);
Z=point_list(:,3);
l=direction(1);
m=direction(2);
n=direction(3);

% notice l, m, n is known
matrix=[-2*point_list,ones(point_number,1)];
U=-sum(point_list.^2,2)+sum(l^2*X.^2 +2*l*m*X.*Y+2*l*n*X.*Z+m^2*Y.^2+2*m*n*Y.*Z +n^2*Z.^2,2)/C;
coeff=(matrix\U); % x0, y0, z0, D

base_point=coeff(1:3)';
x0=base_point(1);
y0=base_point(2);
z0=base_point(3);

D=sum(base_point.^2)+sum(-2*l^2*X*x0+l^2*x0^2-2*l*m*x0*Y+2*l*m*x0*y0-2*l*n*x0*Z+2*l*n*x0*z0-2*l*m*X*y0-2*m^2*Y*y0+m^2*y0^2-2*m*n*y0*Z+2*m*n*y0*z0-2*l*n*X*z0-2*m*n*Y*z0-2*n^2*Z*z0+n^2*z0^2)/C-coeff(4);
D=abs(D); % here have some problem, if not all point are on the cylinder
radius=sqrt(D);

% calculate error
error_list=(sum((point_list-base_point).^2,2)-radius*radius)-...
    sum((direction.*(point_list-base_point)).^2,2)/C;
error=sum(error_list.^2)/point_number;
end
