function S=airEntropy(H,P)
% base pressure and enthalpy calculate entropy
%
% input:
% H (J/kg), P (Pa)
%
% output:
% S (J/(kg*K))
%
% reference: [1] 张志成. 高超声速气动热和热防护[M]. 北京：国防工业出版社, 2003.
%

% [1] Page: 90
P_1=101325;
H_0=4.1868e3; % J/kg
if H <= 1161.8e3 % 167.5 KJ/kg < h <= 1161.8 KJ/kg
    S=3.07311-0.28711*log(P/P_1)+1.00479*log(H/H_0);
else % 1161.8 kJ/kg < h < 62383 kJ/kg
    S=(5.98796-0.195649*log(P/P_1))*exp((H/418.68e3)^0.275/3.548);
end

end
