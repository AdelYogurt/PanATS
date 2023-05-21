function mu=airViscosityTemperature(T)
% base temperature calculate viscosity
%
mu=1.716*1e-5*(T/273.15).^1.5*(273.15+110.4)./(T+110.4);
end