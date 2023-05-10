# D-BAAM
Dual-adsorbent Batch Adsorber Analogue Model (D-BAAM)

This repository contains MATLAB code for executing the D-BAAM process modelling framework presented in [1]. Briefly, the model represents separation of a two-component gas phase mixture by a four-step vacuum swing adsorption process with light product pressurization. The process cycle is intended to collect the heavy component at high purity during the evacuation step. The model code calculates four key performance indicators of the process performance; the purity (Pu) and recovery (Re) of the extracted heavy component, the volumetric working capacity of the system (WC), and the specific energy usage (ET) per tonne of the heavy component which is collected.

The main model function is "DBAAM.m", which calls the D-BAAM model to evaluate the process KPIs. An example of how to use "DBAAM.m" is provided in the file "D-BAAM_example.m". The modelling presented herin in adapted from that of [2]. Full details of the model formulation and performance calculations is provided in [1].

We also provide a database containing the isotherm parameters of CO2/N2 on 76 previously studied materials for post-combustion CO2 capture. The isotherm of each material is given in the dual-site Langmuir (DSL) form. The adsorption isotherms presented in the database were digitsied from [3].

[1] Ward et al. (2023): "Assessment of dual-adsorbent beds for CO2 capture by equilibrium-based process design". Separation and Purification Technology. https://doi.org/10.1016/j.seppur.2023.123990

[2] Balashanker et al. (2019): "Analysis of a Batch Adsorber Analogue for Rapid Screening of Adsorbents for Postcombustion CO2 Capture". Industrial & Engineering Chemistry Research (58). https://doi.org/10.1021/acs.iecr.8b05420

[3] Khurana et al. (2016): "Adsorbent Screening for postcombustion CO2 Capture: A Method Relating Equilibrium Isotherm Characteristics to an Optimum Vacuum Swing Adsorption Process Performance". Industrial & Engineering Chemistry Research (55). https://doi.org/10.1021/acs.iecr.5b04531
