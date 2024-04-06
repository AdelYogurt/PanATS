function mu=airViscosity(H,P)
% base enthalpy calculate viscosity
%
mu_1=1.697e-5;
H_0=262.25e3; % J/kg
A_0=114.55e3; % J/kg

% Sutherland equation
mu=mu_1*(H/H_0).^1.5.*(H_0+A_0)./(H+A_0);

% if H < 711.7e3
%     mu=0.296109e-6*H^0.7259;
% elseif H < 1881.1e3
%     mu=1.05311e-6*H^0.5325;
% else
%     b=1.75+3.99551e-2*log(P/101325);
%     mu=58.4182e-6*(0.31469*log(H/4.1868e3)-0.9221)^b;
% end

% mu=1.716*1e-5*(T/273.15).^1.5*(273.15+110.4)./(T+110.4);
end
