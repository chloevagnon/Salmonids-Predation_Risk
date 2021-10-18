# Salmonids-Predation_Risk
Data, functions and R script are provided for applying analyses of the article:
"The vulnerability of the European whitefish (*Coregonus lavaretus*) to the invasive European catfish (*Silurus glanis*) in a large peri-Alpine lake"  

The script **cv_Rcodes_Vulnerability_Salmonidae** is devided into two parts: (i) Computation of the size-based vulnerability and (ii) Computation of the predation risk.  
  
For processing the size-based vulnerability:  
  - **cv_Functions_aNM**: functions from Vagnon *et al.* (2021) for computing the prey body size range of a consumer, infering trophic links, weighting trophic links
  - **Param_reginvert/Param_regvert**: requested parameters for computing the prey body size range
  - **DATA_BodySize**: data of *S. glanis* and *C. lavaretus* body sizes  
    
For processing the predation risk:  
  - **cv_Functions_Vul_Salmo**: functions to find the best parmeter combination for the LogNormal probability density function (this step can take a long time so it is advised to use the parameters already fitted for the article for a quick application), to compute the depth-matching and the basal metabolism 
  - **Param**: parameters computed from the functions Param_Depth_Matching in cv_Functions_Vul_Salmo.R
  - **Depth_SIL**: data of *S. glanis* depth
  - **Depth_C_0**: data of 0+ whitefish depths
  - **Depth_SIL**: data of the other whitefish lifestages depths
  - **TempMonth**: data of temperatures for computing metabolism and the final predation risk metric


For any supplementary questions, contact us to chloe.vagnon@gmail.com
