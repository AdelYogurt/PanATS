function [max_heat_flow]=solveModelHypersonicHeat()
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

dimension=user_model.dimension;

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
gama=user_model.GAMA_VALUE;
Re=user_model.REYNOLDS_NUMBER;
T_w=user_model.MARKER_ISOTHERMAL;

% load data from inviscid and streamline result
delta_list=output_inviscid.delta_list;
P_list=output_inviscid.P_list;
streamline_len_list=output_streamline.streamline_len_list;
attachment_list=output_streamline.attachment_list;

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
H_0=(airEnthalpy(T_1)+V_1_sq/2); % Free-stream enthalpy J/kg equal to H_s
H_w=airEnthalpy(T_w); % wall air enthalpy J/kg
H_D=(33867*0.78+15320*0.22)*1e3; % Mean dissociation enthalpy of air J/kg
Re_x_tri=10^(5.37+0.2326*Ma_1-0.004015*Ma_1*Ma_1); % transition Reynolds number

% Sutherland method
miu_1=airViscosity(airEnthalpy(T_1));
miu_w=airViscosity(airEnthalpy(T_w));

% initialize result sort array
rou_e_list=cell(length(marker_list),1);
V_e_list=cell(length(marker_list),1);
miu_e_list=cell(length(marker_list),1);
H_e_list=cell(length(marker_list),1);
P_2_list=cell(length(marker_list),1);
rou_2_list=cell(length(marker_list),1);
T_2_list=cell(length(marker_list),1);

Re_x_list=cell(length(marker_list),1);
Q_list=cell(length(marker_list),1);

% calculate boundary layer air paremeter and local Reynolds number
for monitor_index=1:length(MARKER_MONITORING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITORING(monitor_index),marker_list);

    rou_e_list_marker=zeros(marker_list(marker_index).element_number,1);
    V_e_list_marker=zeros(marker_list(marker_index).element_number,1);
    miu_e_list_marker=zeros(marker_list(marker_index).element_number,1);
    H_e_list_marker=zeros(marker_list(marker_index).element_number,1);
    P_2_list_marker=zeros(marker_list(marker_index).element_number,1);
    rou_2_list_marker=zeros(marker_list(marker_index).element_number,1);
    T_2_list_marker=zeros(marker_list(marker_index).element_number,1);

    Re_x_list_marker=zeros(marker_list(marker_index).element_number,1);
    Q_list_marker=zeros(marker_list(marker_index).element_number,1);

    for element_index=1:marker_list(marker_index).element_number
        element=marker_element(element_index);
        point_index_list=element.point_index_list;
        streamline_len=streamline_len_list{marker_index}(element_index);

        % boundary layer air paremeter
        P_e=P_list{marker_index}(element_index); % surface pressure equal to boundary layer pressure
        theta=delta_list{marker_index}(element_index); % actually is attack angle

        [rou_e,V_e,miu_e,H_e,P_2,rou_2,T_2,H_2]=aerodynamicBoundaryLayer...
            (rou_1,V_1,T_1,P_1,Ma_1,gama,P_e,H_0,theta,...
            gama_sub,gama_plus);

        % local Reynolds number
        Re_x=rou_e*V_e*streamline_len/miu_e;

        rou_e_list_marker(element_index)=rou_e;
        V_e_list_marker(element_index)=V_e;
        miu_e_list_marker(element_index)=miu_e;
        H_e_list_marker(element_index)=H_e;
        P_2_list_marker(element_index)=P_2;
        rou_2_list_marker(element_index)=rou_2;
        T_2_list_marker(element_index)=T_2;

        Re_x_list_marker(element_index)=Re_x;
    end

    rou_e_list{marker_index}=rou_e_list_marker;
    V_e_list{marker_index}=V_e_list_marker;
    miu_e_list{marker_index}=miu_e_list_marker;
    H_e_list{marker_index}=H_e_list_marker;
    P_2_list{marker_index}=P_2_list_marker;
    rou_2_list{marker_index}=rou_2_list_marker;
    T_2_list{marker_index}=T_2_list_marker;

    Re_x_list{marker_index}=Re_x_list_marker;
    Q_list{marker_index}=Q_list_marker;
end
max_heat_flow=0;

% calculate attachment element by Detra-Kemp-Riddell method
for attachment_index=1:length(attachment_list)
    marker_index=attachment_list(attachment_index,1);
    element_index=attachment_list(attachment_index,2);
    element=marker_list(marker_index).element_list(element_index);

    % stagnation point air parameter parpare
    P_es=P_list{marker_index}(element_index); % stagnation point pressure

    rou_es=rou_e_list{marker_index}(element_index);
    V_es=V_e_list{marker_index}(element_index);
    miu_es=miu_e_list{marker_index}(element_index);
    H_es=H_e_list{marker_index}(element_index);
    P_2s=P_2_list{marker_index}(element_index);
    rou_2s=rou_2_list{marker_index}(element_index);
    T_2s=T_2_list{marker_index}(element_index);

    P_ws=P_es;
    rou_ws=airDensity(H_w,P_ws);

%     % load nearby point coordination to fit ellipsoid
%     point_arou_index_list=[];
%     point_index_list=element.point_index_list;
%     for point_index=1:length(point_index_list)
%         vertex_index=point_index_list(point_index);
%         edge=edge_list(vertex_index);
%         point_arou_index_list=[point_arou_index_list;edge.vertex_ref_list(1:edge.edge_number)];
%     end
%     point_arou_index_list=unique(point_arou_index_list);
%     point_arou_list=point_list(point_arou_index_list,1:3);
% 
%     % fitting sphere
%     [~,radius_s_list]=calFitEllipsoid(point_arou_list);
%     radius_s=mean(radius_s_list);radius_s_min=min(radius_s_list);radius_s_max=max(radius_s_list);
    
%     Q_s=calTauber(radius_s_min);
%     Q_s=calKempRiddell(radius_s_min);
%     Q_s=calDetraKempRiddell(radius_s_min,radius_s_min);

%     % stagnation parameter
%     T_0=T_1*(gama_sub/2*Ma_1*Ma_1+1);
%     P_0=P_1*(gama_sub/2*Ma_1*Ma_1+1)^(gama/gama_sub);
%     rou_0=rou_1*(gama_sub/2*Ma_1*Ma_1+1)^(1/gama_sub);
%     miu_0=airViscosity(airEnthalpy(T_0));
%     T_es=T_1*(gama_sub/2*Ma_1*Ma_1+1);
%     rou_es=rou_1*(gama_plus*Ma_1*Ma_1)/(gama_sub*Ma_1*Ma_1+2);
%     P_es=P_1*((gama_plus^gama_plus*Ma_1^(2*gama))/...
%         (2^gama_plus*gama*Ma_1*Ma_1-gama_sub))^(1/gama_sub);
%     P_w=P_1*(2*gama*Ma_1*Ma_1-gama_sub)/gama_plus;
%     rou_w=rou_es;
%     miu_es=airViscosity(airEnthalpy(T_es));

%     Q_s=calFayRiddell(radius_s_min);

    % load around element
    V_e_arou_list=[];
    center_point_arou_list=[];
    point_index_list=element.point_index_list;
    point_number=length(point_index_list);
    for point_index=1:point_number
        vertex_index=point_index_list(point_index);
        if point_index == point_number
            vertex_ref_index=point_index_list(1);
        else
            vertex_ref_index=point_index_list(point_index+1);
        end

        edge_arou=edge_list(vertex_ref_index);
        index=edge_arou.getRefIndex(vertex_index);

        if ~isempty(index)
            element_arou=edge_arou.element_ref_list(index);

            V_e_arou_list=[V_e_list{element_arou.marker_index}(element_arou.element_index);V_e_arou_list];
            center_point_arou_list=[element_arou.center_point;center_point_arou_list];
        end
    end

    % calculate du_e__ds by around V_e
    dist_list=sqrt(sum((center_point_arou_list-element.center_point).^2,2));
    du_e__ds=sum(V_e_arou_list./dist_list)/point_number;
    Q_s=0.763*Pr^-0.6*((rou_ws*miu_w)/(rou_es*miu_es))^0.1*...
        sqrt(rou_es*miu_es*du_e__ds)*(1+(Le^a-1)*H_D/H_es)*(H_es-H_w);

    if ~isreal(Q_s)
        disp('?');
    end

    Re_x_list{marker_index}(element_index)=0;
    Q_list{marker_index}(element_index)=Q_s;

    if max_heat_flow < Q_s
        max_heat_flow=Q_s;
    end

end

% calculate heat flow by plate reference enthalpy method
for monitor_index=1:length(MARKER_MONITORING)
    [marker_element,marker_index]=getMarkerElement(MARKER_MONITORING(monitor_index),marker_list);

    for element_index=1:marker_list(marker_index).element_number
        % if element is cross by attachment line, do not need to calculate
        if marker_element(element_index).attachment
            continue;
        end

        % load data
        rou_e=rou_e_list{marker_index}(element_index);
        V_e=V_e_list{marker_index}(element_index);
        miu_e=miu_e_list{marker_index}(element_index);
        H_e=H_e_list{marker_index}(element_index);
        P_2=P_2_list{marker_index}(element_index);
        rou_2=rou_2_list{marker_index}(element_index);
        T_2=T_2_list{marker_index}(element_index);

        Re_x=Re_x_list{marker_index}(element_index);

        % reference air paremeter
        H_r=H_e+r_l*V_e*V_e/2;
        H_ref=0.19*H_r+0.23*H_e+0.58*H_w;
        miu_ref=airViscosity(H_ref);
        rou_ref=airDensity(H_ref,P_e);

        if Re_x <= 1e5
            % laminar flow
            Q=0.332*Pr^(-3/2)*rou_e*V_e*Re_x^-0.5*(H_r-H_w)*...
                (rou_ref*miu_ref/rou_e/miu_e)^0.5;
        elseif Re_x < 1e7
            % transition flow
            Q=0.0296*Pr^(-3/2)*rou_e*V_e*Re_x^-0.2*(H_r-H_w)*...
                (rou_ref/rou_e)^0.8*(miu_ref/miu_e)^0.2;
        else
            % turbulent flow
            Q=0.185*Pr^(-3/2)*rou_ref*V_e*(log(Re_x)/log(10))^2.584*(H_r-H_w)*...
                (rou_ref/rou_e)^0.8*(miu_ref/miu_e)^0.2;
        end

        if ~isreal(Q)
            disp('?');
        end

        Q_list{marker_index}(element_index)=Q;
    end
end

output_heat.rou_e_list=rou_e_list;
output_heat.V_e_list=V_e_list;
output_heat.miu_e_list=miu_e_list;
output_heat.H_e_list=H_e_list;
output_heat.P_2_list=P_2_list;
output_heat.rou_2_list=rou_2_list;
output_heat.T_2_list=T_2_list;

output_heat.Re_x_list=Re_x_list;
output_heat.Q_list=Q_list;

user_model.output_heat=output_heat;


    function Q_s=calFayRiddell(radius_s)
        % Fay-Riddell method
        du_e__ds=sqrt(2*(P_es-P_1)/rou_es)/radius_s;
        Q_s=0.763*Pr^-0.6*((rou_ws*miu_w)/(rou_es*miu_es))^0.1*...
            sqrt(rou_es*miu_es*du_e__ds)*(1+(Le^a-1)*H_D/H_es)*(H_es-H_w);
    end

    function Q_s=calTauber(radius_s)
        % Tauber method
        Q_s=1.83e-4/sqrt(radius_s)*sqrt(rou_1)*V_1^3*(1-H_w/H_es);
    end

    function Q_s=calKempRiddell(radius_s)
        % Kemp-Riddell method
        % Q_s=1.103317e8/sqrt(radius_s)*...
        %     (rou_1/rou_sl)^0.5*(V_1/7925)^3.15*(H_es-H_w)/(H_es-3.0145e5);
        % modify
        Q_s=1.103317e8/sqrt(radius_s)*sqrt(2)*...
            (rou_1/rou_sl)^0.5*(V_1/7925)^3.5*(H_es-H_w)/(H_es-3.0145e5);
    end

    function Q_s=calDetraKempRiddell(radius_s_min,radius_s_max)
        % Detra-Kemp-Riddell method
        Q_s=1.1037e8/sqrt(radius_s_min)*sqrt(1.1+0.9*sqrt(radius_s_min/radius_s_max))*...
            (rou_1/rou_sl)^0.5*(V_1/7925)^3.5*(H_es-H_w)/(H_es-3.0145e5);
    end
end

function [center_point,radius]=calFitEllipsoid(point_list)

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

function [rou_e,V_e,miu_e,H_e,P_2,rou_2,T_2,H_2]=aerodynamicBoundaryLayer...
    (rou_1,V_1,T_1,P_1,Ma_1,gama,P_e,H_0,theta,...
    gama_sub,gama_plus)
% function to calculate outer boundary of boundary layer parameter
% P_e is pressure of outer boundary of boundary layer, equal to pressure of surface
% theta is surface attack angle
% P_2, rou_2, T_2 is the air parameter after shock wave
%
[P_2__P_1,rou_2__rou_1,T_2__T_1]=aerodynamicShockWaveSimple...
    (gama,Ma_1,theta,...
    gama_sub,gama_plus);

P_2=P_1*P_2__P_1;
rou_2=rou_1*rou_2__rou_1;
T_2=T_1*T_2__T_1;
H_2=gama/gama_sub*P_2/rou_2;

% base on isentropic relationship obtain density
rou_e=(P_e/P_2)^(1/gama)*rou_2; % notic from after shock to boundary out is isentropic

% % base on ideal gas equation obtain enthalpy, gama/gama_sub=cp
H_e=gama/gama_sub*P_e/rou_e;
% H_e=(0.213833*P_e/101325/rou_e)^(1/0.972)*1755.52e3;
% if H_e > 1755.52e3
%     a=0.718+1.38974e-2*log(P_e/101325);
%     H_e=(0.213833*P_e/101325/rou_e)^(1/a)*1755.52e3;
% end

% base on energy equation obtain velocity
V_e=sqrt(2*(H_0-H_e));

% base on sutherland to calculate viscosity coefficient
miu_e=airViscosity(H_e);
    function [P_2__P_1,rou_2__rou_1,T_2__T_1]=aerodynamicShockWaveSimple...
            (gama,Ma_1,theta,...
            gama_sub,gama_plus)
        % function to calculate parameter after shock wave
        %
        if abs(theta-pi/2) < 1e-12
            % normal shock
            beta=pi/2;
        else
            % oblique shock wave
            beta=asin(sqrt(functionBetaTheta(gama,Ma_1,theta)));
            beta=beta(2);
            if ~isreal(beta)
                beta=real(beta);
            end
        end
        
        sin_beta_sq=(sin(beta))^2;
        Ma_1_sq=Ma_1*Ma_1;
        
        % calculate parameter
        if abs(beta-pi/2) < 1e-12
            % normal shock
            P_2__P_1=2*gama/gama_plus*Ma_1_sq-gama_sub/gama_plus;
            rou_2__rou_1=gama_plus*Ma_1_sq/(gama_sub*Ma_1_sq+2);
            T_2__T_1=(gama_sub/gama_plus)^2*(2*gama*Ma_1_sq/gama_sub-1)*...
                (2/gama_sub/Ma_1_sq+1);
        else
            % oblique shock wave
            P_2__P_1=2*gama/gama_plus*Ma_1_sq*sin_beta_sq-gama_sub/gama_plus;
            rou_2__rou_1=gama_plus*Ma_1_sq*sin_beta_sq/(gama_sub*Ma_1_sq*sin_beta_sq+2);
            T_2__T_1=(gama_sub/gama_plus)^2*(2*gama*Ma_1_sq*sin_beta_sq/gama_sub-1)*...
                (2/gama_sub/Ma_1_sq/sin_beta_sq+1);
        end
        function sin_beta_sq=functionBetaTheta(gama,Ma_1,theta)
            % function to get sin(beta)^2 by theta
            %
            tan_theta_sq__=tan(theta)^2;
            Ma_1_sq__=Ma_1*Ma_1;
            Ma_1_qu__=Ma_1_sq__*Ma_1_sq__;
            gama_plus__=gama+1;
            c0=1;
            c1=-(tan_theta_sq__*(Ma_1_qu__*gama_plus__*gama_plus__/4+Ma_1_sq__*gama_plus__+1)+...
                (2*Ma_1_sq__+1));
            c2=(tan_theta_sq__*(Ma_1_qu__*gama_plus__+2*Ma_1_sq__)+(Ma_1_qu__+2*Ma_1_sq__));
            c3=-(tan_theta_sq__*Ma_1_qu__+Ma_1_qu__);
            sin_beta_sq=roots([c3 c2 c1 c0]);
            sin_beta_sq=sin_beta_sq(1:2);
        end
        function tan_theta=functionThetaBeta(gama,Ma_1,beta)
            % function to get tan(theta) by beta
            %
            sin_beta__=sin(beta);
            sin_beta_sq__=sin_beta__*sin_beta__;
            tan_beta__=tan(beta);
            Ma_1_sq__=Ma_1*Ma_1;
            tan_theta=(Ma_1_sq__*sin_beta_sq__-1)/...
                (Ma_1_sq__*((gama+1)/2-sin_beta_sq__)+1)/tan_beta__;
        end
    end
end

function [beta,theta,Ma_2,...
    P_2__P_1,rou_2__rou_1,T_2__T_1,P_12__P_11]=aerodynamicShockWave...
    (gama,Ma_1,beta,theta)
% function to calculate aerodynamic parameter after oblique shock wave
% beta is oblique shock wave angle, theta is pointed wedge angle
% gama is fluid specific heat ratio, Ma_1 is fluid mach number
% p01, p02 is total pressure
%
if nargin < 4
    theta=[];
    if nargin < 3
        error('aerodynamicShockWave: lack input angle');
    end
end

if Ma_1 < 1
    error('aerodynamicShockWave: Ma_1 less than 1');
end

if isempty(beta)
    % input theta obtain beta
    if abs(theta-pi/2) < 1e-12
        % normal shock
        beta=pi/2;
    else
        % oblique shock wave
        beta=asin(sqrt(functionBetaTheta(gama,Ma_1,theta)));
        beta=beta(2);
    end
end

if isempty(theta)
    % input beta obtain theta
    if abs(beta-pi/2) < 1e-12
        theta=pi/2;
    else
        tan_theta=functionThetaBeta(gama,Ma_1,beta);
        theta=atan(tan_theta);
    end
end

gama_sub=gama-1;
gama_plus=gama+1;
sin_beta_sq=(sin(beta))^2;
Ma_1_sq=Ma_1*Ma_1;

% calculate parameter
if abs(beta-pi/2) < 1e-12
    % normal shock
    P_2__P_1=2*gama/gama_plus*Ma_1_sq-gama_sub/gama_plus;
    rou_2__rou_1=gama_plus*Ma_1_sq/(gama_sub*Ma_1_sq+2);
    T_2__T_1=(gama_sub/gama_plus)^2*(2*gama*Ma_1_sq/gama_sub-1)*...
        (2/gama_sub/Ma_1_sq+1);
    P_12__P_11=(2*gama*Ma_1_sq/gama_plus-gama_sub/gama_plus)^(-1/gama_sub)*...
        (gama_plus*Ma_1_sq/(gama_sub*Ma_1_sq+2))^(gama/gama_sub);
    Ma_2_sq=(Ma_1_sq+2/gama_sub)/(2*gama*Ma_1_sq/gama_sub-1);
    Ma_2=sqrt(Ma_2_sq);
else
    % oblique shock wave
    P_2__P_1=2*gama/gama_plus*Ma_1_sq*sin_beta_sq-gama_sub/gama_plus;
    rou_2__rou_1=gama_plus*Ma_1_sq*sin_beta_sq/(gama_sub*Ma_1_sq*sin_beta_sq+2);
    T_2__T_1=(gama_sub/gama_plus)^2*(2*gama*Ma_1_sq*sin_beta_sq/gama_sub-1)*...
        (2/gama_sub/Ma_1_sq/sin_beta_sq+1);
    P_12__P_11=(2*gama*Ma_1_sq*sin_beta_sq/gama_plus-gama_sub/gama_plus)^(-1/gama_sub)*...
        (gama_plus*Ma_1_sq*sin_beta_sq/(gama_sub*Ma_1_sq*sin_beta_sq+2))^(gama/gama_sub);
    Ma_2_sq=(Ma_1_sq+2/gama_sub)/(2*gama*Ma_1_sq*sin_beta_sq/gama_sub-1)+...
        2/gama_sub*Ma_1_sq*cos(beta)^2/(Ma_1_sq*sin_beta_sq+2/gama_sub);
    Ma_2=sqrt(Ma_2_sq);
end

    function sin_beta_sq=functionBetaTheta(gama,Ma_1,theta)
        % function to get sin(beta)^2 by theta
        %
        tan_theta_sq__=tan(theta)^2;
        Ma_1_sq__=Ma_1*Ma_1;
        Ma_1_qu__=Ma_1_sq__*Ma_1_sq__;
        gama_plus__=gama+1;
        c0=1;
        c1=-(tan_theta_sq__*(Ma_1_qu__*gama_plus__*gama_plus__/4+Ma_1_sq__*gama_plus__+1)+...
            (2*Ma_1_sq__+1));
        c2=(tan_theta_sq__*(Ma_1_qu__*gama_plus__+2*Ma_1_sq__)+(Ma_1_qu__+2*Ma_1_sq__));
        c3=-(tan_theta_sq__*Ma_1_qu__+Ma_1_qu__);
        sin_beta_sq=roots([c3 c2 c1 c0]);
        sin_beta_sq=sin_beta_sq(1:2);
    end
    function tan_theta=functionThetaBeta(gama,Ma_1,beta)
        % function to get tan(theta) by beta
        %
        sin_beta__=sin(beta);
        sin_beta_sq__=sin_beta__*sin_beta__;
        tan_beta__=tan(beta);
        Ma_1_sq__=Ma_1*Ma_1;
        tan_theta=(Ma_1_sq__*sin_beta_sq__-1)/...
            (Ma_1_sq__*((gama+1)/2-sin_beta_sq__)+1)/tan_beta__;
    end
end

function miu=airViscosity(H,P)
% base enthalpy calculate viscosity
%
miu_1=1.697e-5;
H_0=262.25e3; % J/kg
A_0=114.55e3; % J/kg

% Sutherland equation
miu=miu_1*(H/H_0).^1.5.*(H_0+A_0)./(H+A_0);

% if H < 711.7e3
%     miu=0.296109e-6*H^0.7259;
% elseif H < 1881.1e3
%     miu=1.05311e-6*H^0.5325;
% else
%     b=1.75+3.99551e-2*log(P/101325);
%     miu=58.4182e-6*(0.31469*log(H/4.1868e3)-0.9221)^b;
% end

% miu=1.716*1e-5*(T/273.15).^1.5*(273.15+110.4)./(T+110.4);
end

function S=airEntropy(H,P)
% base pressure and enthalpy calculate entropy
% entropy is kJ/(kg*K), H is kJ/kg
%
P_1=101325;
H_0=4.1868e3; % J/kg
if H <= 1161.8e3 % 167.5kJ/kg<h<1161.8kJ/kg
    S=3.07311-0.28711*log(P/P_1)+1.00479*log(H/H_0);
else % 1161.8kJ/kg<h<62383kJ/kg
    S=(5.98796-0.195649*log(P/P_1))*exp((H/418.68e3)^0.275/3.548);
end
end

function rou=airDensity(H,P)
% base pressure and enthalpy calculate entropy
% H is kJ/kg
P_1=101325;
if H <= 1755.52e3 % 167.5kJ/kg<h<=1755.52kJ/kg
    a=0.972;
else % 1755.52kJ/kg<h<35026kJ/kg
    a=0.718+1.38974e-2*log(P/P_1);
end
rou=0.213833*(P/P_1)*(H/1755.52e3)^-a;
end

function H=airEnthalpy(T)
% base tempareture calculate enthalpy
%
% W=28.97;
% R=8314.472;
% H=R/W*(3.51715*T+2.5041*1e-4*T^2+5.5076*1e-7*T^3-1.7197*1e-10*T^4);

% 200K - 1200K
% C_list=[1.009432005343459e+03,0.071868378789092,1.580696789340697e-04,-4.935587676354850e-08];
% H=C_list*[T;T.*T;T.*T.*T;T.^4];

% 50K - 3000K
C_list=[1.817160e3,9.890504e2,-9.595592e-3,1.041469e-4,-4.433065e-8,5.879263e-12];
H=C_list*[ones(1,length(T));T;T.*T;T.*T.*T;T.^4;T.^5];
end

