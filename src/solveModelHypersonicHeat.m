function [max_heat_flow]=solveModelHypersonicHeat()
% plate reference enthalpy method to calculate hypersonic aircraft
%
% point_list is coordinate of all node
% element_list contain element(element_type, node_index1, node_index2, ...)
% marker_list contain make{marker_name,marker_element_number,marker_element}
% marker_element contain contain element(element_type, node_index1, node_index2, ...)
%
% note P_e == P_w because of
% pressure of outer boundary of boundary layer is equal to pressure of surface
%
% copyright Adel 2022.11
%
global g_model
    
marker_element_number=size(HATS_element_list,1);
dimension=3;

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

g_geometry.rou_1=rou_1;
g_geometry.V_1=V_1;
g_geometry.T_1=T_1;
g_geometry.P_1=P_1;
g_geometry.Ma_1=Ma_1;
g_geometry.gama=gama;

% solve prepare
H_0=(airEnthalpy(T_1)+V_1_sq/2); % Free-stream enthalpy J/kg equal to H_s
H_w=airEnthalpy(T_w); % wall air enthalpy J/kg
H_d=(33867*0.78+15320*0.22)*1e3; % Mean dissociation enthalpy of air J/kg

streamline_length_list=streamline_output.streamline_length_list;
delta_list=inviscid_output.delta_list;
P_list=inviscid_output.P_list;

% Sutherland method
miu_1=airViscosity(airEnthalpy(T_1));
miu_w=airViscosity(airEnthalpy(T_w));

% stagnation point air parameter parpare
P_es=P_list(1); % stagnation point pressure
P_ws=P_es;
theta=pi/2; % actually is attack angle
[rou_es,V_es,miu_es,H_es,P_2,rou_2,T_2]=aerodynamicBoundaryLayer...
    (rou_1,V_1,T_1,P_1,Ma_1,gama,P_es,H_0,theta,...
    gama_sub,gama_plus);
rou_w_s=airDensity(H_w,P_ws);

% % stagnation parameter
% T_0=T_1*(gama_sub/2*Ma_1*Ma_1+1);
% P_0=P_1*(gama_sub/2*Ma_1*Ma_1+1)^(gama/gama_sub);
% rou_0=rou_1*(gama_sub/2*Ma_1*Ma_1+1)^(1/gama_sub);
% miu_0=airViscosity(airEnthalpy(T_0));
%
% T_s=T_1*(gama_sub/2*Ma_1*Ma_1+1);
% rou_s=rou_1*(gama_plus*Ma_1*Ma_1)/(gama_sub*Ma_1*Ma_1+2);
% P_s=P_1*((gama_plus^gama_plus*Ma_1^(2*gama))/...
%     (2^gama_plus*gama*Ma_1*Ma_1-gama_sub))^(1/gama_sub);
% P_w=P_1*(2*gama*Ma_1*Ma_1-gama_sub)/gama_plus;
% rou_w=rou_s;
% miu_s=airViscosity(airEnthalpy(T_s));

% % Fay-Riddell method
% du_e__ds=sqrt(2*(P_es-P_w)/rou_es)/radius_s;
% Q_s_FR=0.763*Pr^-0.6*(rou_w_s*miu_w/rou_es/miu_es)^0.1*...
%     sqrt(rou_es*miu_es*du_e__ds)*(1+(Le^a-1)*H_d/H_es)*(H_es-H_w);



if nargin < 10
    % calculate stagnation point heat flow density
    root_element=HATS_element_list(1);
    children_element=HATS_element_list(root_element.out_index_list);
    stagnation_point_index_list=[root_element.point_index_list,[children_element.point_index_list]];
    stagnation_point_index_list=unique(stagnation_point_index_list);
    stagnation_point_list=g_Point(stagnation_point_index_list,1:3);
    
    % fitting sphere
    [~,radius_s]=fitSphere(stagnation_point_list);
    
    % Detra-Kemp-Riddell method
    Q_s=1.1037e8/sqrt(radius_s)*sqrt(1.1+0.9)*...
        (rou_1/rou_sl)^0.5*(V_1/7925)^3.5*(H_es-H_w)/(H_es-3.0145e5);
else
    % Detra-Kemp-Riddell method
    Q_s=1.1037e8/sqrt(D_x)*sqrt(1.1+0.9*sqrt(D_x/D_z))*...
        (rou_1/rou_sl)^0.5*(V_1/7925)^3.5*(H_es-H_w)/(H_es-3.0145e5);
end

% % Kemp-Riddell method
% Q_s_KR=1.103317e8/sqrt(radius_s)*...
%     (rou_1/rou_sl)^0.5*(V_1/7900)^3.15*(H_es-H_w)/(H_es-airEnthalpy(300));

max_heat_flow=Q_s;

% calculate heat flow by plate reference enthalpy method
Q_list=ones(marker_element_number,1)*Q_s;
Re_x_list=zeros(marker_element_number,1);
for element_index=2:marker_element_number
    streamline_length=streamline_length_list(element_index);
    
    % boundary layer air paremeter
    P_e=P_list(element_index); % surface pressure equal to boundary layer pressure
    theta=delta_list(element_index); % actually is attack angle
    [rou_e,V_e,miu_e,H_e,P_2,rou_2,T_2]=aerodynamicBoundaryLayer...
        (rou_1,V_1,T_1,P_1,Ma_1,gama,P_e,H_0,theta,...
        gama_sub,gama_plus);
    
    % reference air paremeter 
    H_r=H_e+r_l*V_e*V_e/2;
    H_ref=0.19*H_r+0.23*H_e+0.58*H_w;
    miu_ref=airViscosity(H_ref);
    rou_ref=airDensity(H_ref,P_e);
    
    % Local Reynolds number
    Re_x=rou_e*V_e*streamline_length/miu_e;
    
    % calculate heat flow density
    Q=0.332*Pr^(-3/2)*rou_e*V_e*Re_x^-0.5*(H_r-H_w)*...
        (rou_ref*miu_ref/rou_e/miu_e)^0.5;
    
%     if Re_x < 2540
%         % laminar flow
%         Q=0.332*Pr^(-3/2)*rou_e*V_e*Re_x^-0.5*(H_r-H_w)*...
%             (rou_ref*miu_ref/rou_e/miu_e)^0.5;
%     else
%         % turbulent flow
%         Q=0.185*Pr^(-3/2)*rou_ref*V_e*(H_r-H_w)/...
%             (log(Re_x*(rou_ref*miu_ref/rou_e/miu_e)))^2.584;
%     end
    
    if length(Q) ~= 1
        disp('error');
    end
    
    Q_list(element_index)=min(Q_list(element_index),Q);
    Re_x_list(element_index)=Re_x;
end

heat_output.Q_list=Q_list;
heat_output.Re_x_list=Re_x_list;

    function [center_point,radius]=fitSphere(point_list)
        % fitting sphere base on point list
        % g_Point is point_number x dimension matrix
        %
        [point_number__,dimension__]=size(point_list);
        average_point=sum(point_list(:,1:3),1)/point_number__;
        point_differ_list=point_list-average_point;
        matrix__=[
            sum(point_differ_list(:,1).^2),sum(point_differ_list(:,1).*point_differ_list(:,2)),sum(point_differ_list(:,1).*point_differ_list(:,3));
            sum(point_differ_list(:,2).*point_differ_list(:,1)),sum(point_differ_list(:,2).^2),sum(point_differ_list(:,2).*point_differ_list(:,3));
            sum(point_differ_list(:,3).*point_differ_list(:,1)),sum(point_differ_list(:,3).*point_differ_list(:,2)),sum(point_differ_list(:,3).^2);];
        colume__=[
            sum(point_differ_list(:,1).^3+point_differ_list(:,1).*point_differ_list(:,2).^2+point_differ_list(:,1).*point_differ_list(:,3).^2)/2;
            sum(point_differ_list(:,2).*point_differ_list(:,1).^2+point_differ_list(:,2).^3+point_differ_list(:,2).*point_differ_list(:,3).^2)/2;
            sum(point_differ_list(:,3).*point_differ_list(:,1).^2+point_differ_list(:,3).*point_differ_list(:,2).^2+point_differ_list(:,3).^3)/2;];
        center_point=matrix__\colume__;
        radius=sqrt(sum((point_differ_list(:,1)-center_point(1)).^2+...
            (point_differ_list(:,2)-center_point(2)).^2+...
            (point_differ_list(:,3)-center_point(3)).^2)/point_number__);
        center_point=center_point+average_point';
    end
end

function [rou_e,V_e,miu_e,H_e,P_2,rou_2,T_2]=aerodynamicBoundaryLayer...
    (rou_1,V_1,T_1,P_1,Ma_1,gama,P_e,H_0,theta,...
    gama_sub,gama_plus)
% function to calculate outer boundary of boundary layer parameter
% P_e is pressure of outer boundary of boundary layer, equal to pressure of surface
% theta is surface attack angle
%
[P_2__P_1,rou_2__rou_1,T_2__T_1]=aerodynamicShockWaveSimple...
    (gama,Ma_1,theta,...
    gama_sub,gama_plus);

P_2=P_1*P_2__P_1;
rou_2=rou_1*rou_2__rou_1;
T_2=T_1*T_2__T_1;

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
                beta=pi/2;
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