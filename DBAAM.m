function [Pu_c,Re_c,WC,E_n] = BLAM(Para1,Para2,PH,PINT,PL,T,w,m,V,e,yc_feed)
% Calculate the process KPIs using the BLAM model

W_BLO = 0; % Blwodown work [kWh]
W_EVA = 0; % Evacuation work [kWh]
R = 8.3145; % Universal gas constant [J/mol/K]
k = 1.4; % Adiabatic constant [-]

% Solve the non-linear system of material balance equations for the ADS and LPP steps
% x = [N_LPP, yc_d, N_feed, N_raff]
x0 = [0.2,0.01,1,1]; % Initial guess
options = optimoptions('fsolve', 'Display', 'off'); % Solver options
x = fsolve(@(x) LPPADS(x, yc_feed, Para1, Para2, PH, PL, T, w, m, V, e), x0, options); % Solve system of non-linear equations

% Solve ODE to obtain P(y) for BLO and EVAC steps
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6); % Solver options
[P, yc] = ode45(@(P,yc) ODEfunc(Para1, Para2, P, T, yc, w, m, V, e), (100000:-100:PL-100), 0.15, options); % Solve ODE
idx1 = find( P <= PINT, 1);
yc_b = yc(idx1);
idx2 = find( P <= PL, 1);
yc_c = yc(idx2);

% Rectangle rule for integrating work done by vacuum pump in BLO and EVAC steps
for i = 1:(idx1 - 1)
    [qcs_a1,qns_a1] = DSL(Para1,P(i),T,yc(i));
    [qcs_a2,qns_a2] = DSL(Para2,P(i),T,yc(i));
    [qcs_b1,qns_b1] = DSL(Para1,P(i+1),T,yc(i+1));
    [qcs_b2,qns_b2] = DSL(Para2,P(i+1),T,yc(i+1));
    [qc_a,qn_a] = Mixing(qcs_a1,qcs_a2,qns_a1,qns_a2,w);
    [qc_b,qn_b] = Mixing(qcs_b1,qcs_b2,qns_b1,qns_b2,w);
    Eta = 0.8 * 19.55 * (P(i)/100000) ./ (1 + 19.55 * (P(i)/100000)); % isentropic efficiency
    deltaN = -((P(i + 1) - P(i)) * V * e ./ (R * T) + m * (qc_b + qn_b - qc_a - qn_a));
    Constcoe = (1 / Eta) * (k / (k - 1)) * R * T;
    W_BLO = W_BLO + Constcoe * ((1 ./ (P(i)/100000)) .^ ((k - 1) / k) - 1) * deltaN;
end

for i = idx1:(idx2 - 1)
    [qcs_a1,qns_a1] = DSL(Para1,P(i),T,yc(i));
    [qcs_a2,qns_a2] = DSL(Para2,P(i),T,yc(i));
    [qcs_b1,qns_b1] = DSL(Para1,P(i+1),T,yc(i+1));
    [qcs_b2,qns_b2] = DSL(Para2,P(i+1),T,yc(i+1));
    [qc_a,qn_a] = Mixing(qcs_a1,qcs_a2,qns_a1,qns_a2,w);
    [qc_b,qn_b] = Mixing(qcs_b1,qcs_b2,qns_b1,qns_b2,w);
    Eta = 0.8 * 19.55 * (P(i)/100000) ./ (1 + 19.55 * (P(i)/100000));
    deltaN = -((P(i + 1) - P(i)) * V * e ./ (R * T) + m * (qc_b + qn_b - qc_a - qn_a));
    Constcoe = (1 / Eta) * (k / (k - 1)) * R * T;
    W_EVA = W_EVA + Constcoe * ((1 ./ (P(i)/100000)) .^ ((k - 1) / k) - 1) * deltaN;
end

% Loading at stage beta for adsorbent 1 and adsorbent 2
[qcs_b1, qns_b1] = DSL(Para1, PINT, T, yc_b);
[qcs_b2, qns_b2] = DSL(Para2, PINT, T, yc_b);

% Loading at stage gamma for adsorbent 1 and adsorbent 2
[qcs_c1, qns_c1] = DSL(Para1, PL, T, yc_c);
[qcs_c2, qns_c2] = DSL(Para2, PL, T, yc_c);

% Loading at stage beta and gamma for layered bed
[qc_b, qn_b] = Mixing(qcs_b1, qcs_b2, qns_b1, qns_b2, w);
[qc_c, qn_c] = Mixing(qcs_c1, qcs_c2, qns_c1, qns_c2, w);

% Calculate KPIs
M_c = 44 / (1e6); % Molecular weight of CO2 [tonne/mol]
CO2_cap = (yc_b * PINT - yc_c * PL) * V * e ./ (R * T) + m * (qc_b - qc_c); % Amount of CO2 captured [mol]

Pu_c = 100*((yc_b * PINT - yc_c * PL) * V * e ./ (R * T) + m * (qc_b - qc_c)) ./ ((PINT - PL) * V * e ./ (R * T) + m * (qc_b + qn_b - qc_c - qn_c)); % Purity [%]
Re_c = 100*((yc_b * PINT - yc_c * PL) * V * e ./ (R * T) + m * (qc_b - qc_c)) ./ (x(3) * yc_feed); % Recovery [%]
WC = ((yc_b * PINT - yc_c * PL) * V * e ./ (R * T) + m * (qc_b - qc_c)) ./ (V * (1 -e)); % Working capacity [mol/m3]
E_n = (W_BLO + W_EVA) / (CO2_cap * M_c) * 2.777778e-7; % Energy usage [kWh/tonne]

end