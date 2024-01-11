function [beta,theta,Ma_2,...
    T_2__T_1,P_2__P_1,rho_2__rho_1,P_02__P_01]=aerodynamicShockWave...
    (gamma,Ma_1,beta,theta)
% function to calculate aerodynamic parameter after oblique shock wave
% beta is oblique shock wave angle, theta is pointed wedge angle
% gamma is fluid specific heat ratio, Ma_1 is fluid mach number
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
        beta=asin(sqrt(functionBetaTheta(gamma,Ma_1,theta)));
        beta=beta(2);
    end
end

if isempty(theta)
    % input beta obtain theta
    if abs(beta-pi/2) < 1e-12
        theta=pi/2;
    else
        tan_theta=functionThetaBeta(gamma,Ma_1,beta);
        theta=atan(tan_theta);
    end
end

gamma_sub=gamma-1;
gamma_plus=gamma+1;
sin_beta_sq=(sin(beta))^2;
Ma_1_sq=Ma_1*Ma_1;

% calculate parameter
if abs(beta-pi/2) < 1e-12
    % normal shock
    T_2__T_1=(gamma_sub/gamma_plus)^2*(2*gamma*Ma_1_sq/gamma_sub-1)*...
        (2/gamma_sub/Ma_1_sq+1);
    P_2__P_1=2*gamma/gamma_plus*Ma_1_sq-gamma_sub/gamma_plus;
    rho_2__rho_1=gamma_plus*Ma_1_sq/(gamma_sub*Ma_1_sq+2);
    P_02__P_01=(2*gamma*Ma_1_sq/gamma_plus-gamma_sub/gamma_plus)^(-1/gamma_sub)*...
        (gamma_plus*Ma_1_sq/(gamma_sub*Ma_1_sq+2))^(gamma/gamma_sub);
    Ma_2_sq=(Ma_1_sq+2/gamma_sub)/(2*gamma*Ma_1_sq/gamma_sub-1);
    Ma_2=sqrt(Ma_2_sq);
else
    % oblique shock wave
    T_2__T_1=(gamma_sub/gamma_plus)^2*(2*gamma*Ma_1_sq*sin_beta_sq/gamma_sub-1)*...
        (2/gamma_sub/Ma_1_sq/sin_beta_sq+1);
    P_2__P_1=2*gamma/gamma_plus*Ma_1_sq*sin_beta_sq-gamma_sub/gamma_plus;
    rho_2__rho_1=gamma_plus*Ma_1_sq*sin_beta_sq/(gamma_sub*Ma_1_sq*sin_beta_sq+2);
    P_02__P_01=(2*gamma*Ma_1_sq*sin_beta_sq/gamma_plus-gamma_sub/gamma_plus)^(-1/gamma_sub)*...
        (gamma_plus*Ma_1_sq*sin_beta_sq/(gamma_sub*Ma_1_sq*sin_beta_sq+2))^(gamma/gamma_sub);
    Ma_2_sq=(Ma_1_sq+2/gamma_sub)/(2*gamma*Ma_1_sq*sin_beta_sq/gamma_sub-1)+...
        2/gamma_sub*Ma_1_sq*cos(beta)^2/(Ma_1_sq*sin_beta_sq+2/gamma_sub);
    Ma_2=sqrt(Ma_2_sq);
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
