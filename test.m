clc;
clear;
close all hidden;

load('matlab.mat')

theta_list=0:0.01:pi/2;
rou_e_list=zeros(1,length(theta_list));
for theta_index=1:length(theta_list)
    theta=theta_list(theta_index);
    rou_e_list(theta_index)=aerodynamicBoundaryLayer...
        (rou_1,V_1,T_1,P_1,Ma_1,gama,P_e,H_0,theta,...
        gama_sub,gama_plus);
end
line(theta_list,rou_e_list);

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


% radius_s=0.0095;
% 
% du_e__ds=sqrt(2*(P_es-P_1)/rou_es)/radius_s;
% Q_s=0.763*Pr^-0.6*(rou_es*miu_e)^0.4*(rou_ws*miu_w)^0.1*(1+(Le^a-1)*H_D/H_es)*...
%     sqrt(du_e__ds)*(H_es-H_w)
% 
% Q_s=1.103317e8/sqrt(radius_s)*sqrt(2)*...
%             (rou_1/rou_sl)^0.5*(V_1/7925)^3.5*(H_es-H_w)/(H_es-3.0145e5)
% 
% rou_es*miu_e
% rou_miu=calRouMiu(P_es,H_es)
% rou_ws*miu_w
% rou_miu=calRouMiu(P_ws,H_w)
% 
% function rou_miu=calRouMiu(P,H)
% rou_b=515.3788184*2.498e-3;
% miu_b=47.88026*3.584e-7;
% P_b=47.8802589*2117;
% H_b=0.0929*2.119e8;
% 
% rou_miu=(0.225*rou_b*miu_b*(P/P_b)^0.992)/...
%     (1-1.0213*(1-(H/H_b)^0.3329));
% 
% end

