function H=airEnthalpy(T)
% base tempareture calculate enthalpy
%
% input:
% T (K)
%
% output:
% H (J/kg)
%
% reference: [1] Miller C G, Iii. Comparison of thin-film resistance
% heat-transfer gages with thin-skin transient calorimeter gages in
% conventional hypersonic wind tunnels[R]. NASA-TM-83197, 1981.
%
%

% [1] Page: 17
% W=28.97;
% R=8314.472;
% H=R/W*(3.51715*T+2.5041*1e-4*T^2+5.5076*1e-7*T^3-1.7197*1e-10*T^4);

% [1] Page: 17
T=T(:);
% 200K - 1200K
C_list=[1.009432005343459e+03,0.071868378789092,1.580696789340697e-04,-4.935587676354850e-08];
H=C_list*[T;T.*T;T.*T.*T;T.^4];

% 
% 50K - 3000K
% C_list=[1.817160e3,9.890504e2,-9.595592e-3,1.041469e-4,-4.433065e-8,5.879263e-12];
% H=C_list*[ones(1,length(T));T;T.*T;T.*T.*T;T.^4;T.^5];
end