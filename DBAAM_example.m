%%
% A routine demonstrating how to call the BLAM model to calculate process KPIs
% of purity, recovery, working capacity and energy usage

clear all; close all; clc;



%%
%======================================
%========== INPUT PARAMETERS ==========
%======================================

yc_feed = 0.15;   % Mole fraction of CO2 in feed [-]
PH = 1e+5;        % High pressure [Pa]
PINT = 0.15e+5;   % Intermediate pressure [Pa]
PL = 0.03e+5;     % Low pressure [Pa]
T = 298.15;       % Temperature [K]
w = 0.5;          % Loading ratio of adsorbent 1 [-]
m = 1;            % Mass of adsorbent [kg]
e = 0.37;         % Void fraction [-]
rho_1 = 1130;     % Particle density of adsorbent 1 [kg/m3]
rho_2 = 799.5;    % Particle density of adsorbent 2 [kg/m3]
V = m*w/(1-e)/rho_1 + m*(1-w)/(1-e)/rho_2; % Volume of system [m3]

% Isotherm parameters (extended DSL model)
% [b0_CO2, d0_CO2, dUb_CO2, dUd_CO2, qsb_CO2, qsd_CO2, b0_N2, d0_N2, dUb_N2, dUd_N2, qsb_N2, qsd_N2]
% b0 or d0 in [m3/mol]
% dUb or dUd in [J/mol]
% qsb or qsd in [mol/kg]
Para1 = {8.65e-7 2.63e-8 -36600 -35700 3.09 2.54 2.5e-6  0 -15800 0 5.84 0}; % DSL parameters of adsorbent 1 (zeolite 13X)
Para2 = {9.4e-06 1.04e-05 -25610 -17550 0.59 7.51 1.81e-03 1.72e-12 -8670 -44900 0.16 41.3}; % DSL parameters of adsorbent 2 (Activated carbon)



%%
%===================================
%========== EVALUATE KPIs ==========
%===================================

[Pu_c, Re_c, WC, E_n] = DBAAM(Para1, Para2, PH, PINT, PL, T, w, m, V, e, yc_feed)

% Pu_c: purity of CO2 [%]
% Re_c: recovery of CO2 [%]
% WC: Working capacity [mol/m3]
% E_n: specific energy usage [kWh/tonne]