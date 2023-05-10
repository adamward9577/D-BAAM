function [dyc] = ODEfunc(Para1, Para2, P, T, yc, w, m, V, e)
% ODE RHS gradient function

R = 8.3145; % Universal gas constant [J/mol/K]

[bco1, dco1, Ubc1, Udc1, qsbc1, qsdc1, bno1, dno1, Ubn1, Udn1, qsbn1, qsdn1] = deal(Para1{:});
[bco2, dco2, Ubc2, Udc2, qsbc2, qsdc2, bno2, dno2, Ubn2, Udn2, qsbn2, qsdn2] = deal(Para2{:});

bc1 = bco1 * exp(-Ubc1 / (R * T));
dc1 = dco1 * exp(-Udc1 / (R * T));
bn1 = bno1 * exp(-Ubn1 / (R * T));
dn1 = dno1 * exp(-Udn1 / (R * T));
bc2 = bco2 * exp(-Ubc2 / (R * T));
dc2 = dco2 * exp(-Udc2 / (R * T));
bn2 = bno2 * exp(-Ubn2 / (R * T));
dn2 = dno2 * exp(-Udn2 / (R * T));

dqcs1_dP = qsbc1 * bc1 * yc * R * T ./ ((bc1 * yc - bn1 * yc + bn1) * P + R * T) .^ 2 + qsdc1 * dc1 * yc * R * T ./ ((dc1 * yc - dn1 * yc + dn1) * P + R * T) .^ 2;
dqcs2_dP = qsbc2 * bc2 * yc * R * T ./ ((bc2 * yc - bn2 * yc + bn2) * P + R * T) .^ 2 + qsdc2 * dc2 * yc * R * T ./ ((dc2 * yc - dn2 * yc + dn2) * P + R * T) .^ 2;

dqns1_dP = qsbn1 * bn1 * (1 - yc) * R * T ./ ((bc1 * yc - bn1 * yc + bn1) * P + R * T) .^ 2 + qsdn1 * dn1 * (1 - yc) * R * T ./ ((dc1 * yc - dn1 * yc + dn1) * P + R * T) .^ 2;
dqns2_dP = qsbn2 * bn2 * (1 - yc) * R * T ./ ((bc2 * yc - bn2 * yc + bn2) * P + R * T) .^ 2 + qsdn2 * dn2 * (1 - yc) * R * T ./ ((dc2 * yc - dn2 * yc + dn2) * P + R * T) .^ 2;

dqcs1_dyc = qsbc1 * bc1 * P * (bn1 * P + R * T) ./ ((bc1 * P - bn1 * P) * yc + bn1 * P + R * T) .^ 2 + qsdc1 * dc1 * P * (dn1 * P + R * T) ./ ((dc1 * P - dn1 * P) * yc + dn1 * P + R * T) .^ 2;
dqcs2_dyc = qsbc2 * bc2 * P * (bn2 * P + R * T) ./ ((bc2 * P - bn2 * P) * yc + bn2 * P + R * T) .^ 2 + qsdc2 * dc2 * P * (dn2 * P + R * T) ./ ((dc2 * P - dn2 * P) * yc + dn2 * P + R * T) .^ 2;

dqns1_dyc = - qsbn1 * bn1 * P * (bc1 * P + R * T) ./ ((bc1 * P - bn1 * P) * yc + bn1 * P + R * T) .^ 2 - qsdn1 * dn1 * P * (dc1 * P + R * T) ./ ((dc1 * P - dn1 * P) * yc + dn1 * P + R * T) .^ 2;
dqns2_dyc = - qsbn2 * bn2 * P * (bc2 * P + R * T) ./ ((bc2 * P - bn2 * P) * yc + bn2 * P + R * T) .^ 2 - qsdn2 * dn2 * P * (dc2 * P + R * T) ./ ((dc2 * P - dn2 * P) * yc + dn2 * P + R * T) .^ 2;

dqc_dP = dqcs1_dP * w + dqcs2_dP * (1 - w);
dqn_dP = dqns1_dP * w + dqns2_dP * (1 - w);
dqc_dyc = dqcs1_dyc * w + dqcs2_dyc * (1 - w);
dqn_dyc = dqns1_dyc * w + dqns2_dyc * (1 - w);

a1 = V * e ./ (R * T) + m * (dqc_dP + dqn_dP);
a2 = yc * V * e ./ (R * T) + m * dqc_dP;
f1 = m * (dqc_dyc + dqn_dyc);
f2 = P * V * e ./ (R * T) + m * dqc_dyc;

dyc = (a1 * yc - a2) ./ (f2 - f1 * yc);

end