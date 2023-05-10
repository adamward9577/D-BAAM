function [qc, qn] = Mixing(qcs1, qcs2, qns1, qns2, w)
% Layered bed loading

qc = qcs1 * w + qcs2 * (1 - w); % Loading of CO2 [mol/kg]
qn = qns1 * w + qns2 * (1 - w); % Loading of N2 [mol/kg]

end