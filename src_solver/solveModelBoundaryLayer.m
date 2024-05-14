function solveModelBoundaryLayer()
% solve boundary layer parameter
%
% note P_e == P_w because of
% pressure of outer boundary of boundary layer is equal to pressure of surface
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

% heat calculate need inviscid and streamline result
output_inviscid=user_model.output_inviscid;
output_streamline=user_model.output_streamline;

T_inf=config.FREESTREAM_TEMPERATURE;
P_inf=config.FREESTREAM_PRESSURE;
Ma_inf=config.MACH_NUMBER;
gamma=config.GAMMA_VALUE;
% Re=config.REYNOLDS_NUMBER;
T_w=config.MARKER_ISOTHERMAL;

% load data from inviscid and streamline result
theta_list=output_inviscid.theta_list;
P_list=output_inviscid.P_list;

streamline_len_list=output_streamline.streamline_len_list;

R=287.0955;
Pr=0.71;
gamma_plus=gamma+1;
gamma_sub=gamma-1;

% free flow parameter
rho_inf=P_inf/R/T_inf;
a_inf=sqrt(gamma*R*T_inf);
V_inf=a_inf*Ma_inf;
V_inf_sq=V_inf*V_inf;
q_inf=rho_inf*V_inf_sq/2;

% solve prepare
Re_x_tri=10^(5.37+0.2326*Ma_inf-0.004015*Ma_inf*Ma_inf); % transition Reynolds number

% temperature recovery coefficient
r_l_lam=Pr^(1/2); % laminar flow
r_l_turb=Pr^(1/3); % turbulence flow

% solve prepare
H_0=(airEnthalpy(T_inf)+V_inf_sq/2); % Free-stream enthalpy J/kg equal to H_s
H_w=airEnthalpy(T_w); % wall air enthalpy J/kg

elem_num=length(element_list);
% initialize result sort array
T_p_list=zeros(elem_num,1);
P_p_list=zeros(elem_num,1);
rho_p_list=zeros(elem_num,1);

rho_e_list=zeros(elem_num,1);
V_e_list=zeros(elem_num,1);
mu_e_list=zeros(elem_num,1);
H_e_list=zeros(elem_num,1);

rho_ref_list=zeros(elem_num,1);
mu_ref_list=zeros(elem_num,1);
H_ref_list=zeros(elem_num,1);

H_r_list=zeros(elem_num,1);

Re_x_list=zeros(elem_num,1);
Re_x_ref_list=zeros(elem_num,1);

elem_num=length(element_list);
% calculate boundary layer air paremeter and
% local Reynolds number, local reference Reynolds number
for elem_idx=1:elem_num
    elem=element_list(elem_idx);
    streamline_len=streamline_len_list(elem_idx);

    if streamline_len == 0
        % for stagnation element or attachmenet line element
        % use average streamline_len in element
        Node_idx=elem.Node_idx;
        node_num=elem.node_num;
        center_point=sum(point_list(Node_idx,1:3),1)/node_num;
        point_len=sqrt(sum((point_list(Node_idx,1:3)-center_point).^2,2));
        streamline_len=sum(point_len)/node_num;
    end

    % boundary layer air paremeter
    P_e=P_list(elem_idx); % surface pressure equal to boundary layer pressure
    theta=theta_list(elem_idx); % actually is attack angle

    [T_p,P_p,rho_p,H_p,rho_e,V_e,mu_e,H_e]=aerodynamicBoundaryLayer...
        (Ma_inf,T_inf,P_inf,rho_inf,V_inf,gamma,P_e,H_0,theta,...
        gamma_sub,gamma_plus);

    if V_e == 0, V_e=V_inf;end

    % local Reynolds number
    Re_x=rho_e*V_e*streamline_len/mu_e;

    % reference air paremeter
    if Re_x < Re_x_tri
        % laminar flow
        H_r=H_e+r_l_lam*V_e*V_e/2;
    else
        % turbulent flow
        H_r=H_e+r_l_turb*V_e*V_e/2;
    end
    H_ref=0.19*H_r+0.23*H_e+0.58*H_w;
    mu_ref=airViscosity(H_ref);
    rho_ref=airDensity(H_ref,P_e);

    % local reference Reynolds number
    Re_x_ref=Re_x*(rho_ref*mu_e)/(rho_e*mu_ref);

    T_p_list(elem_idx)=T_p;
    P_p_list(elem_idx)=P_p;
    rho_p_list(elem_idx)=rho_p;

    rho_e_list(elem_idx)=rho_e;
    V_e_list(elem_idx)=V_e;
    mu_e_list(elem_idx)=mu_e;
    H_e_list(elem_idx)=H_e;

    rho_ref_list(elem_idx)=rho_ref;
    mu_ref_list(elem_idx)=mu_ref;
    H_ref_list(elem_idx)=H_ref;

    H_r_list(elem_idx)=H_r;

    Re_x_list(elem_idx)=Re_x;
    Re_x_ref_list(elem_idx)=Re_x_ref;
end

output_boulay.T_p_list=T_p_list;
output_boulay.P_p_list=P_p_list;
output_boulay.rho_p_list=rho_p_list;

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

if config.INFORMATION
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
H_e=gamma/gamma_sub*P_e^(1-1/gamma)*P_2^(1/gamma)/rho_2; % H_e=gamma/gamma_sub*P_e/rho_e;

% base on energy equation obtain velocity
V_e=sqrt(2*(H_0-H_e));

if H_0 < H_e, V_e=0;end

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
                beta=real(beta);
                % beta=pi/2;
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

