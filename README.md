# Salmonids-Predation_Risk
Data, functions and R scipt are provided for applying analyses of the article:
"The vulnerability of whitefish (Coregonus lavaretus) to the invasive European catfish (Silurus glanis) in a large peri-Alpine lake"  

The script "cv_Rcodes_Vulnerability_Salmonidae" is devided into two parts: (i) Computation of the size-based vulnerability and (ii) Computation of the predation risk.  
  
For processing the size-based vulnerability:  
  - cv_Functions_aNM: functions for computing the prey body size range of a consumer, infering trophic links, weighting trophic links
  - Param_reginvert/Param_regvert: paramters necessary for computing the prey body size range
  - DATA_BodySize: data of S. glanis and C. lavaretus body sizes  
    
For processing the predation risk:  
  - cv_Functions_Vul_Salmo: functions to find the best parmeter combination for the LogNormal probability density function and to compute the depth-matching 
  - Param_LogNormal: parameters computed from the functions in cv_Functions_Vul_Salmo
  - Depth_SIL: data of S. glanis depth
  - Depth_C_0: data of 0+ whitefish depths
  - Depth_SIL: data of the other whitefish lifestages depths
