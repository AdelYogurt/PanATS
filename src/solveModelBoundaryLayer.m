function solveModelBoundaryLayer()
% solve boundary layer parameter
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

% load data from inviscid and streamline result
theta_list=output_inviscid.theta_list;
P_list=output_inviscid.P_list;
streamline_len_list=output_streamline.streamline_len_list;

R=287.0955;
Pr=0.71;
r_l=Pr^(1/3); % temperature recovery coefficient
gamma_plus=gamma+1;
gamma_sub=gamma-1;

% free flow parameter
rho_1=P_1/R/T_1;
a_1=sqrt(gamma*R*T_1);
V_1=a_1*Ma_1;
V_1_sq=V_1*V_1;
q_1=rho_1*V_1_sq/2;

% solve prepare
H_0=(airEnthalpy(T_1)+V_1_sq/2); % Free-stream enthalpy J/kg equal to H_s
H_w=airEnthalpy(T_w); % wall air enthalpy J/kg

% initialize result sort array
T_2_list=cell(length(marker_list),1);
P_2_list=cell(length(marker_list),1);
rho_2_list=cell(length(marker_list),1);

rho_e_list=cell(length(marker_list),1);
V_e_list=cell(length(marker_list),1);
mu_e_list=cell(length(marker_list),1);
H_e_list=cell(length(marker_list),1);

rho_ref_list=cell(length(marker_list),1);
mu_ref_list=cell(length(marker_list),1);
H_ref_list=cell(length(marker_list),1);

H_r_list=cell(length(marker_list),1);

Re_x_list=cell(length(marker_list),1);
Re_x_ref_list=cell(length(marker_list),1);

% calculate boundary layer air paremeter and 
% local Reynolds number, local reference Reynolds number
for monitor_index=1:length(MARKER_MONITORING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITORING(monitor_index),marker_list);

    T_2_list_marker=zeros(marker_list(marker_index).element_number,1);
    P_2_list_marker=zeros(marker_list(marker_index).element_number,1);
    rho_2_list_marker=zeros(marker_list(marker_index).element_number,1);
    
    rho_e_list_marker=zeros(marker_list(marker_index).element_number,1);
    V_e_list_marker=zeros(marker_list(marker_index).element_number,1);
    mu_e_list_marker=zeros(marker_list(marker_index).element_number,1);
    H_e_list_marker=zeros(marker_list(marker_index).element_number,1);

    rho_ref_list_marker=zeros(marker_list(marker_index).element_number,1);
    mu_ref_list_marker=zeros(marker_list(marker_index).element_number,1);
    H_ref_list_marker=zeros(marker_list(marker_index).element_number,1);

    H_r_list_marker=zeros(marker_list(marker_index).element_number,1);

    Re_x_list_marker=zeros(marker_list(marker_index).element_number,1);
    Re_x_ref_list_marker=zeros(marker_list(marker_index).element_number,1);

    for element_index=1:marker_list(marker_index).element_number
        element=marker_element(element_index);
        streamline_len=streamline_len_list{marker_index}(element_index);

        if streamline_len == 0
            % for stagnation element or attachmenet line element
            % use average streamline_len in element
            point_index_list=element.point_index_list;
            center_point=sum(point_list(point_index_list,1:3),1)/length(point_index_list);
            point_len=sqrt(sum(point_list(point_index_list,1:3)-center_point,2));
            streamline_len=sum(point_len)/length(point_index_list);
        end

        % boundary layer air paremeter
        P_e=P_list{marker_index}(element_index); % surface pressure equal to boundary layer pressure
        theta=theta_list{marker_index}(element_index); % actually is attack angle

        [T_2,P_2,rho_2,H_2,rho_e,V_e,mu_e,H_e]=aerodynamicBoundaryLayer...
            (Ma_1,T_1,P_1,rho_1,V_1,gamma,P_e,H_0,theta,...
            gamma_sub,gamma_plus);

%         V_e=sqrt(sum((V_1*element.surface_flow).^2));

        % local Reynolds number
        Re_x=rho_e*V_e*streamline_len/mu_e;

        % reference air paremeter
        H_r=H_e+r_l*V_e*V_e/2;
        H_ref=0.19*H_r+0.23*H_e+0.58*H_w;
        mu_ref=airViscosity(H_ref);
        rho_ref=airDensity(H_ref,P_e);

        % local reference Reynolds number
        Re_x_ref=Re_x*(rho_ref*mu_e)/(rho_e*mu_ref);

        T_2_list_marker(element_index)=T_2;
        P_2_list_marker(element_index)=P_2;
        rho_2_list_marker(element_index)=rho_2;

        rho_e_list_marker(element_index)=rho_e;
        V_e_list_marker(element_index)=V_e;
        mu_e_list_marker(element_index)=mu_e;
        H_e_list_marker(element_index)=H_e;

        rho_ref_list_marker(element_index)=rho_ref;
        mu_ref_list_marker(element_index)=mu_ref;
        H_ref_list_marker(element_index)=H_ref;

        H_r_list_marker(element_index)=H_r;

        Re_x_list_marker(element_index)=Re_x;
        Re_x_ref_list_marker(element_index)=Re_x_ref;
    end

    T_2_list{marker_index}=T_2_list_marker;
    P_2_list{marker_index}=P_2_list_marker;
    rho_2_list{marker_index}=rho_2_list_marker;

    rho_e_list{marker_index}=rho_e_list_marker;
    V_e_list{marker_index}=V_e_list_marker;
    mu_e_list{marker_index}=mu_e_list_marker;
    H_e_list{marker_index}=H_e_list_marker;

    rho_ref_list{marker_index}=rho_ref_list_marker;
    mu_ref_list{marker_index}=mu_ref_list_marker;
    H_ref_list{marker_index}=H_ref_list_marker;

    H_r_list{marker_index}=H_r_list_marker;

    Re_x_list{marker_index}=Re_x_list_marker;
    Re_x_ref_list{marker_index}=Re_x_ref_list_marker;
end

output_boulay.T_2_list=T_2_list;
output_boulay.P_2_list=P_2_list;
output_boulay.rho_2_list=rho_2_list;

output_boulay.rho_e_list=rho_e_list;
output_boulay.V_e_list=V_e_list;
output_boulay.mu_e_list=mu_e_list;
output_boulay.H_e_list=H_e_list;

output_boulay.rho_ref_list=rho_ref_list;
output_boulay.mu_ref_list=mu_ref_list;
output_boulay.H_ref_list=H_ref_list;

output_boulay.H_r_list=H_r_list;

output_boulay.Re_x_list=Re_x_list;
output_boulay.Re_x_ref_list=Re_x_ref_list;

user_model.output_boulay=output_boulay;

if user_model.INFORMATION
    fprintf('solveModelBoundaryLayer: boundary layer parameter solve done!\n');
end

end

function [T_2,P_2,rho_2,H_2,rho_e,V_e,mu_e,H_e]=aerodynamicBoundaryLayer...
    (Ma_1,T_1,P_1,rho_1,V_1,gamma,P_e,H_0,theta,...
    gamma_sub,gamma_plus)
% function to calculate outer boundary of boundary layer parameter
% P_e is pressure of outer boundary of boundary layer, equal to pressure of surface
% theta is surface attack angle
% P_2, rho_2, T_2 is the air parameter after shock wave
%
[T_2__T_1,P_2__P_1,rho_2__rho_1]=aerodynamicShockWaveSimple...
    (gamma,Ma_1,theta,...
    gamma_sub,gamma_plus);

P_2=P_1*P_2__P_1;
rho_2=rho_1*rho_2__rho_1;
T_2=T_1*T_2__T_1;
H_2=gamma/gamma_sub*P_2/rho_2;

% base on isentropic relationship obtain density
rho_e=(P_e/P_2)^(1/gamma)*rho_2; % notic from after shock to boundary out is isentropic

% % base on ideal gas equation obtain enthalpy, gamma/gamma_sub*R=cp
H_e=gamma/gamma_sub*P_e/rho_e;
% H_e=(0.213833*P_e/101325/rho_e)^(1/0.972)*1755.52e3;
% if H_e > 1755.52e3
%     a=0.718+1.38974e-2*log(P_e/101325);
%     H_e=(0.213833*P_e/101325/rho_e)^(1/a)*1755.52e3;
% end

% base on energy equation obtain velocity
V_e=sqrt(2*(H_0-H_e));

if H_0 < H_e
    V_e=0;
end

if ~isreal(V_e)
    disp('?');
end

% base on sutherland to calculate viscosity coefficient
mu_e=airViscosity(H_e);
    function [T_2__T_1,P_2__P_1,rho_2__rho_1]=aerodynamicShockWaveSimple...
            (gamma,Ma_1,theta,...
            gamma_sub,gamma_plus)
        % function to calculate parameter after shock wave
        %
        if abs(theta-pi/2) < 1e-12
            % normal shock
            beta=pi/2;
        else
            % oblique shock wave
            beta=asin(sqrt(functionBetaTheta(gamma,Ma_1,theta)));
            beta=beta(2);
            if ~isreal(beta)
%                 beta=real(beta);
                beta=pi/2;
            end
        end
        
        sin_beta_sq=(sin(beta))^2;
        Ma_1_sq=Ma_1*Ma_1;
        
        % calculate parameter
        if abs(beta-pi/2) < 1e-12
            % normal shock
            P_2__P_1=2*gamma/gamma_plus*Ma_1_sq-gamma_sub/gamma_plus;
            rho_2__rho_1=gamma_plus*Ma_1_sq/(gamma_sub*Ma_1_sq+2);
            T_2__T_1=(gamma_sub/gamma_plus)^2*(2*gamma*Ma_1_sq/gamma_sub-1)*...
                (2/gamma_sub/Ma_1_sq+1);
        else
            % oblique shock wave
            P_2__P_1=2*gamma/gamma_plus*Ma_1_sq*sin_beta_sq-gamma_sub/gamma_plus;
            rho_2__rho_1=gamma_plus*Ma_1_sq*sin_beta_sq/(gamma_sub*Ma_1_sq*sin_beta_sq+2);
            T_2__T_1=(gamma_sub/gamma_plus)^2*(2*gamma*Ma_1_sq*sin_beta_sq/gamma_sub-1)*...
                (2/gamma_sub/Ma_1_sq/sin_beta_sq+1);
        end

        function sin_beta_sq=functionBetaTheta(gamma,Ma_1,theta)
            % function to get sin(beta)^2 by theta
            %
            tan_theta_sq__=tan(theta)^2;
            Ma_1_sq__=Ma_1*Ma_1;
            Ma_1_qu__=Ma_1_sq__*Ma_1_sq__;
            gama_plus__=gamma+1;
            c0=1;
            c1=-(tan_theta_sq__*(Ma_1_qu__*gama_plus__*gama_plus__/4+Ma_1_sq__*gama_plus__+1)+...
                (2*Ma_1_sq__+1));
            c2=(tan_theta_sq__*(Ma_1_qu__*gama_plus__+2*Ma_1_sq__)+(Ma_1_qu__+2*Ma_1_sq__));
            c3=-(tan_theta_sq__*Ma_1_qu__+Ma_1_qu__);
            sin_beta_sq=roots([c3 c2 c1 c0]);
            sin_beta_sq=sin_beta_sq(1:2);
        end

        function tan_theta=functionThetaBeta(gamma,Ma_1,beta)
            % function to get tan(theta) by beta
            %
            sin_beta__=sin(beta);
            sin_beta_sq__=sin_beta__*sin_beta__;
            tan_beta__=tan(beta);
            Ma_1_sq__=Ma_1*Ma_1;
            tan_theta=(Ma_1_sq__*sin_beta_sq__-1)/...
                (Ma_1_sq__*((gamma+1)/2-sin_beta_sq__)+1)/tan_beta__;
        end

    end

end

