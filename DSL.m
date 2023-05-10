function [qcs, qns] = DSL(Para, P, T, yc)
% Dual-site Langmuir isotherm model

R = 8.3145; % Universal gas constant [J/mol/K]
yn = 1 - yc; % Mole fraction of N2 [-]

[bco, dco, Ubc, Udc, qsbc, qsdc, bno, dno, Ubn, Udn, qsbn, qsdn] = deal(Para{:});

bc = bco * exp(-Ubc / (R * T)); % Equilibrium constant of CO2 at site 1 [m3/mol]
dc = dco * exp(-Udc / (R * T)); % Equilibrium constant of CO2 at site 2 [m3/mol]
bn = bno * exp(-Ubn / (R * T)); % Equilibrium constant of N2 at site 1 [m3/mol]
dn = dno * exp(-Udn / (R * T)); % Equilibrium constant of N2 at site 2 [m3/mol]

Cc = yc * P / (R * T); % Concentration of CO2 [m3/mol]
Cn = yn * P / (R * T); % Concentration of N2 [m3/mol]

qcs = qsbc * bc * Cc ./ (1 + bc * Cc + bn * Cn) + qsdc * dc * Cc ./ (1 + dc * Cc + dn * Cn); % Amount of CO2 adsorbed [mol/kg]
qns = qsbn * bn * Cn ./ (1 + bc * Cc + bn * Cn) + qsdn * dn * Cn ./ (1 + dc * Cc + dn * Cn); % Amount of N2 adsorbed [mol/kg]

end