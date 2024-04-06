function rho=airDensity(H,P)
% base pressure and enthalpy calculate entropy
% H is kJ/kg
P_1=101325;
if H <= 1755.52e3 % 167.5kJ/kg<h<=1755.52kJ/kg
    a=0.972;
else % 1755.52kJ/kg<h<35026kJ/kg
    a=0.718+1.38974e-2*log(P/P_1);
end
rho=0.213833*(P/P_1)*(H/1755.52e3)^-a;
end
