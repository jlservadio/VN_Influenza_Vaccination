# VN_Influenza_Vaccination

Code and data accompanying "Influenza vaccination allocation in tropical settings under constrained resources," available on medRxiv (https://www.medrxiv.org/content/10.1101/2024.02.08.24302551v1) and submitted for peer-reviewed publication. 

## Contents:

**betas_for_start_Full.Rdata** contains 50 time-varying parameter values for the transmission parameter used in the model for this study

**condenseContact.R** contains code to transform demographic data contained in various age groupings into the 10-year age bands used in this study. These include the contact matrices from Prem et al. (2021) as well as influenza hospitalization and fatality ratios from the US CDC, 

**Implement_Vaccine_Strategies_10pct.R** contains code to implement the model and test vaccine allocations under a supply available for 10% of the population. 

**Make_Allocations_LHS.R** contains code for Latin hypercube sampling to generate vaccine allocations for each vaccine supply. 

**model_code_multiplegroup.txt** contains Rcpp code outlining the mathematical model used in this study.

**Params_for_start_full.Rdata** contains all parameterizations of the model considered in this study. Analyses were repeated for 50 parameterizations for robustness. 

**VN_Demographics.Rdata** contains age demographic data for Vietnam.

**VN_US_pops.Rdata** contains age demographic data for Vietnam, the United States, and intermediate demographics used for sensitivity analysis. 
