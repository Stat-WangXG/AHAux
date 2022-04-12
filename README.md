# AHAux

Supplement to ‘Auxiliary summary information synthesis method based on additive hazard model’.

R code and data for analysis in Auxiliary summary information synthesis method based on additive hazard model.

Here is a general explanation of what is contained in the folders.

R

AH.R: Fitting of the additive hazard model based on the estimation procedure proposed by Lin and Ying (1994).
AH_Aux_GMM.R: Fitting of the additive hazard model with auxiliary covariate effects information based on GMM.
AH_EE.R: The details of estimating equation based on the additive hazard model at the individual level. 
Calibration.R: A calibration method for parameter estimation.
gmm_weight.R: The estimation of the weight matrix in the objective function of generalized moments method.
J_test.R: The Sargan-Hansen J-test for the conformity condition.

Data

simdata.RData: A simulated right-censored dataset (25% censoring rate) from the additive hazard model with four covariates.

