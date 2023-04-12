function [T,P,rho,a,miu,g]=getAtmosphereEnvironment(Z)
% base on altitude calculate atmosphere parameter
% Z is altitude m
% return temperature pressure density speed_of_sound acceleration_of_gravity
%
Z=Z/1e3;

R_earth=    6.356766e3;     % earth ratio km
T_SL=       2.8815e2;       % sea level temperature
P_SL=       1.01325e5;      % sea level pressure
rou_SL=     1.2250;         % sea level density
g_SL=       9.80665;        % sea level acceleration of gravity
a_SL=       3.40294e2;      % sea level speed of sound

% geopotential height H
H=Z/(1+Z/R_earth);

% calculate parameter
if Z <= 11.0191
    W=1-H/44.3308;
    T=W*T_SL;
    P=W^5.2553*P_SL;
    rho=W^4.2559*rou_SL;
    
elseif Z <= 20.0631
    W=exp((14.9647-H)/6.3416);
    T=216.650;
    P=1.1953*1e-1*W*P_SL;
    rho=1.5898*1e-1*W*rou_SL;
    
elseif Z <= 32.1619
    W=1+(H-24.9021)/221.552;
    T=221.552*W;
    P=2.5158*1e-2*W^-34.1629*P_SL;
    rho=3.2722*1e-2*W^-35.1829*rou_SL;
    
elseif Z <= 47.3501
    W=1+(H-39.4799)/89.4107;
    T=250.350*W;
    P=2.8338*1e-3*W^-12.2011*P_SL;
    rho=3.2617*1e-3*W^-13.2011*rou_SL;
    
elseif Z <= 51.4125
    W=exp((48.6252-H)/7.9223);
    T=270.650;
    P=8.9155*1e-4*W*P_SL;
    rho=9.4920*1e-4*W*rou_SL;
    
elseif Z <=  71.8020
    W=1-(H-59.4390)/88.2218;
    T=247.021*W;
    P=2.1671*1e-4*W^12.2011*P_SL;
    rho=2.5280*1e-4*W^11.2011*rou_SL;
    
elseif Z <= 86.0000
    W=1-(H-78.0303)/100.2950;
    T=200.590*W;
    P=1.2274*1e-5*W^17.0816*P_SL;
    rho=1.7632*1e-5*W^16.0816*rou_SL;
    
elseif Z <= 91.0000
    W=exp((87.2848-H)/5.4700);
    T=186.870;
    P=(2.2730+1.042*1e-3*H)*1e-6*W*P_SL;
    rho=3.6411*1e-6*W*rou_SL;
    
elseif Z <= 150.0000
    T=-0.0040*Z^3+1.5054*Z^2-177.5620*Z+6.8929*1e-3;
    P=2.532*1e6*exp(-0.1829*Z)+0.1403*exp(-0.03698*Z);
    rho=70.22*1e6*exp(-0.1874*Z)+1.734*1e-5*exp(-0.05828*Z);
    
else
    
    
end

a=20.0468*sqrt(T);
g=g_SL/(1+Z/R_earth)^2;

miu=1.716*1e-5*(T/273.15).^1.5*(273.15+110.4)./(T+110.4);

end