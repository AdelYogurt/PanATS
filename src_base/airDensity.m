function rho=airDensity(H,P)
% base pressure and enthalpy calculate entropy
%
% input:
% H (J/kg), P (Pa)
%
% output:
% rho (kg/m3)
%
% reference:
% [1] 张志成. 高超声速气动热和热防护[M]. 北京：国防工业出版社, 2003.
%

% [1] Page: 90
P_1=101325;
if H <= 1755.52e3 % 167.5 kJ/kg < h <= 1755.52 kJ/kg
    a=0.972;
else % 1755.52 kJ/kg < h < 35026 kJ/kg
    a=0.718+1.38974e-2*log(P/P_1);
end
rho=0.213833*(P/P_1)*(H/1755.52e3)^-a;

end
