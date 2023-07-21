function [F] = LPPADS(x, yc_feed, Para1, Para2, PH, PL, T, w, m, V, e)
% System of non-linear equations for the ADS and LPP steps

R = 8.3145; % Universal gas constant [J/mol/K]

% Solve ODE to obtain P(y)
options = odeset('RelTol',1e-6,'AbsTol',1e-6); % Solver options
[P,yc] = ode45(@(P,yc) ODEfunc(Para1, Para2, P, T, yc, w, m, V, e), [100000:-100:PL-100], 0.15, options); % Solve ODEs
idx = find( P <= PL, 1);
yc_c = yc(idx);

% Loading at stage alpha, delta and gamma for adsorbent 1 and adsorbent 2
[qcs_a1, qns_a1] = DSL(Para1, PH, T, yc_feed);
[qcs_a2, qns_a2] = DSL(Para2, PH, T, yc_feed);
[qcs_c1, qns_c1] = DSL(Para1, PL, T, yc_c);
[qcs_c2, qns_c2] = DSL(Para2, PL, T, yc_c);
[qcs_d1, qns_d1] = DSL(Para1, PH, T, x(2));
[qcs_d2, qns_d2] = DSL(Para2, PH, T, x(2));

% Loading at stage alpha, delta and gamma for layered bed
[qc_a, qn_a] = Mixing(qcs_a1, qcs_a2, qns_a1, qns_a2,w);
[qc_c, qn_c] = Mixing(qcs_c1, qcs_c2, qns_c1, qns_c2,w);
[qc_d, qn_d] = Mixing(qcs_d1, qcs_d2, qns_d1, qns_d2,w);

% RHS functions of system of non-linear algebraic equations
F(1) = (PL - PH) * V * e ./ (R * T) + m * (qc_c + qn_c - qc_d - qn_d) + x(1);
F(2) = (yc_c * PL - x(2) * PH) * V * e ./ (R * T) + m * (qc_c - qc_d) + x(1) * x(2);
F(3) = m * (qc_d + qn_d - qc_a - qn_a) + x(3) - x(4);
F(4) = (x(2) - yc_feed) * PH * V * e ./ (R * T) + m * (qc_d - qc_a) + x(3) * yc_feed - x(4) * x(2);

end