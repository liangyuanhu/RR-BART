# RR-BART We propose an inference-based method, called RR-BART, which leverages the likelihood-based Bayesian machine learning technique,
# Bayesian additive regression trees, and uses Rubin's rule to combine the estimates and variances of the variable importance measures on 
# multiply imputed datasets for variable selection in the presence of MAR data. We conduct a representative simulation study to 
# investigate the practical operating characteristics of RR-BART, and compare it with the bootstrap imputation based methods. 
# We further demonstrate the methods via a case study of risk factors for 3-year incidence of metabolic syndrome among middle-aged women 
# using data from the Study of Women's Health Across the Nation (SWAN). The simulation study suggests that even in complex conditions of 
# nonlinearity and nonadditivity with a large percentage of missingness,  RR-BART can reasonably recover both prediction and variable selection 
# performances,  achievable on the fully observed data. RR-BART provides the best performance that the bootstrap imputation based methods can
# achieve with the optimal selection threshold value. In addition, RR-BART demonstrates a substantially stronger ability of detecting discrete 
# predictors. Furthermore, RR-BART offers substantial computational savings. When implemented on the SWAN data, RR-BART adds to the literature by 
# selecting a set of predictors that had been less commonly identified as risk factors but had substantial biological justifications. 
# The accompanying paper Lin et al. (2022) "A flexible approach for variable selection in large-scale healthcare database studies with 
# missing covariate and outcome data" is forthcoming in BMC Medical Research Methodology. 
