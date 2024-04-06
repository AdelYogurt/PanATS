function Cp=airSpecificHeatPressure(T)
% Predict constant pressure specific heat capacity Cp (J/(kg·K)) w.r.t. temperature T (K).
% Perfect gas assumption (T < 2500K) adopted.
% Values at standard atmosphere (atm) (101325kPa) used.
% 
% Reference: [1] HILSENRATH J. Tables of thermal properties of gases:
% Comprising tables of thermodynamic and transport properties of air,
% argon, carbon dioxide, carbon monoxide, hydrogen, nitrogen, oxygen, and
% steam[M]. US Department of Commerce, National Bureau of Standards, 1955.
%
T=T(:);
C_list=[1.0436e+03,-0.3502,9.0412e-04,-5.8095e-07,1.2648e-10];
Cp=C_list*[ones(1,length(T));T;T.*T;T.*T.*T;T.^4];
end