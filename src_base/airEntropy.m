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
