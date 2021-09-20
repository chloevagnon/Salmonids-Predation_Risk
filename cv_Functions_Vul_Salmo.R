
############################################### Function for computing a modal value ########################################
#                                             
# FROSSARD V., VAGNON C./ USMB / 2021
# This function allows to compute the modal value of a number sequence
# Input :
# x =depth sequence
#
#Output : 
# modal value of the depth sequence x

dmode <- function(x) {
  den <- density(x, kernel=c("gaussian"))
  ( den$x[den$y==max(den$y)] )   
}  



############################################## Function for computing the metabolic rate ####################################
#                                              
# FROSSARD V., VAGNON C. / USMB / 2021
# Compute metabolism in KJ/Kg/day adapted from:
# Brown JH, Gillooly JF, Allen AP, et al. (2004) Toward a Metabolic Theory of Ecology.
# Ecology 85:1771-1789. doi: 10.1890/03-9000.
#
# Inputs: 
#  - Body mass of the predator (in grams) 
#  - water temperature (in °C)
# 
# Outputs: 
# - a vector of the metabolism of the predator (in KJ.kg^-1. day^-1)
# /!\ 1 watt is equal to 1 J/s so we multiplied M_B by 3600 to convert second in hour and then, 24 to convert hour in day.

basal_metabolism <- function(bodymass,temperature) {
  Temp_K<-temperature + 273.15 #K
  e=0.63 #eV
  k=8.62*10^-5 #eV.K^-1
  M_B <- (exp(18.47)*(bodymass^(3/4))*exp(-e/(k*Temp_K))) #Watt
  metab <- ((M_B*3600)*24) # J.day^-1
  metab2<-(metab/1000)*(1000/bodymass) # KJ.kg^-1.day^-1
  
  
  return(metab2)
}

             

####################################### Functions for the Spatial vulnerability ################################################# 
#                                      
# 
# VAGNON C. / USMB / May 2021 
# First, the density curves corresponding to the probability of encountering the predator considering its modal depth is fitted.
# The more representative curve is the one resulting from a Log-Normal probability density function.
# The spatial vulnerability is then calculated from the density curve for a given depth.
#
# To do that, two functions were created : 
# 1) "Param_Spatial_Vul" : find the best association of parameters (meanlog and sdlog) for fittig the log normal density curve
#                          for a given modal depth of a predator so that the depth for which the maximum of the density curve 
#                          is found match the modal depth of the predator. The sdlog is fitted so that modal value +-5meters 
#                          include 75 +- 0.5 % of data.
#                                
# 2) "Spatial_Vul" : find the spatial vulnerability of a prey depending on its depth according to the density curve previously fitted.
                       

# FUNCTION 1 : Param_Spatial_Vul
#
# Inputs :
# - DepthAll = sequence from minimal depth to log(maximum depth) of a predator each 0.01meters 
# - DepthMod = modal value of the depth of a predator for each month (can be one or several values)
# - MEAN = meanlog = value of the log(mean) for the simulation of the lognormal density curve
#                    (a vector of a sequence of values to test)
# - SD = sdlog = value of the log(standard deviation) for the simulation of the lognormal density curve
#                    (can be a fixed value or a vector of values SD = seq(0.05,1,0.01))

# Outputs : 
# A list of length 12 (number of months) countaining : 
#    - 12 data.frame with : depth = DepthAll
#                           Prob = density curve fitted for each depth in DepthAll with the best mean for the depthmode of the month
#    - value of the optimal mean and the optimal sd used in the dlnorm to obtain the density curve of vulnerability

Param_Spatial_Vul<-function(DepthAll,DepthMod,MEAN,SD){
  
  Param<-list(Density=vector(mode = "list", length = length(DepthMod)),Param=matrix(NA,nrow=length(DepthMod),ncol=2)) #Initiate param lists
  
  if(length(SD)>1){
    for (k in 1:length(DepthMod)){ #k is the number of months
    
    inter<-matrix(NA,ncol=4)
    
    # Adjust the mean and the SD for obtaining maximum value of the density curves at the depth mode
    for (j in MEAN){
      inter2<-matrix(NA,ncol=4,nrow=length(SD))
      m=0
      
      for(l in SD){
        m=m+1
        uni_mod<-dlnorm(DepthAll, meanlog = j, sdlog = l, log = FALSE) #LogNormal distribution
        
        depth_mod<-DepthAll[which.max(uni_mod)]# Depth for which the density curve reaches the max
        
        inter2[m,1]<-j
        inter2[m,2]<-abs(DepthMod[k]-depth_mod)
        
        
        FN<-ecdf(rlnorm(10^5, meanlog =j,sdlog =l))# rlnorm 
        PercMin<-FN(round(DepthMod[k]+5,digits=1))*100
        PercMax<-FN(round(DepthMod[k]-5,digits=1))*100
        
        
        inter2[m,3]<-l
        inter2[m,4]<-PercMin-PercMax
      }
      inter<-rbind(inter,inter2)
      
    }
    inter<-na.omit(inter)
    
    # Keep optimal parameters 
    op<-inter[which.min(inter[,2]) & inter[,4]<70.5 & inter[,4]>=70,]
    opMEAN<-op[which.min(op[,2]),1]
    opSD<-op[which.min(op[,2]),3]
    
    
    Param[["Density"]][[k]]<-data.frame(depth=DepthAll,Prob=scales::rescale(dlnorm(DepthAll, meanlog = opMEAN, sdlog = opSD, log = FALSE),to=c(0,1)))
    Param[["Param"]][,1][k]<-opMEAN
    Param[["Param"]][,2][k]<-opSD 
    }
  }

  if(length(SD)==1){
    for (k in 1:12){ #k is the number of months
      
      inter<-matrix(NA,ncol=2)
      
      for (j in MEAN){
        inter2<-matrix(NA,ncol=4,nrow=length(SD))
        m=0
        
        for(l in SD){
          m=m+1
          uni_mod<-dlnorm(DepthAll, meanlog = j, sdlog = l, log = FALSE) #LogNormal distribution
          
          depth_mod<-DepthAll[which.max(uni_mod)]# Depth for which the density curve is max
          
          inter2[m,1]<-j
          inter2[m,2]<-abs(DepthMod[k]-depth_mod)
          
        }
        inter<-rbind(inter,inter2)
        
      }
      inter<-na.omit(inter)
      opMEAN<-inter[which.min(inter[,2]),]
      opSD<-SD
      
      
      Param[["Density"]][[k]]<-data.frame(depth=DepthAll,Prob=scales::rescale(dlnorm(DepthAll, meanlog = opMEAN, sdlog = opSD, log = FALSE),to=c(0,1)))
      Param[["Param"]][,1][k]<-opMEAN
      Param[["Param"]][,2][k]<-opSD 
    }
  }
  
  
  return(Param)
}




### FUNCTION 2: Spatial_Vul
# --> Find the vulnerability for given lognormal parameters, Depths occupation of prey and months
# Inputs :
# - Param = Parameters in which the density curves for the month are listed that can be: 
#                - a list resulting from the function "Param_Spatial_Vul_SD"
#                - a data.frame of the same length of Month (1 set of parameters per month)
#                - a vector of 2 parameters that remains the same used in each calculation
# - DepthPrey = depth for which the vulnerability is wanted and can be given as : 
#                - one value 
#                - a sequence of values
# - Month = the month number for which the vulnerability is wanted and can be given as :
#                - one value 
#                - a sequence of values
# - DepthAll = sequence from minimal depth to log(maximum depth) of a predator each 0.01meters 
#             (optional argument if Param is given as a data.frame or as a vector)
#
# Outputs : 
# A data.frame countaining : 
#  - Month : number varying between 1 and 12 according to the input parameters
#  - Depth : Depth of the prey for which the vulnerability was calculated according to the input parameters
#  - Vul : the vulnerability values computed for each depth and month according to the input parameters

Spatial_Vul<-function(Param,DepthPrey,Month,DepthAll){
  
  Spatial_Vulnerability<-data.frame(Month=rep(Month,each=length(DepthPrey)),
                                    Depth=rep(DepthPrey),Vul=NA)
  
  if(class(Param)=="list"){
    for(i in Spatial_Vulnerability$Month){
      for(j in Spatial_Vulnerability$Depth){
        Spatial_Vulnerability$Vul[Spatial_Vulnerability$Month==i &Spatial_Vulnerability$Depth==j]<-Param[["Density"]][[i]]$Prob[which(Param[["Density"]][[i]]$depth==j)]
        }
      }
  }
  
  
  if(class(Param)=="data.frame"){
    for(i in Spatial_Vulnerability$Month){
      for(j in Spatial_Vulnerability$Depth){
        
        Dens<-data.frame(Depth=DepthAll,Prob=scales::rescale(dlnorm(DepthAll,meanlog = Param[i,1],sdlog = Param[i,2],log = FALSE),to=c(0,1)))
        
        Spatial_Vulnerability$Vul[Spatial_Vulnerability$Month==i &Spatial_Vulnerability$Depth==j]<-Dens$Prob[which(Dens$Depth==j)]
      }
    }
  }
  
  
  if(class(Param)=="numeric"){
    for(i in Spatial_Vulnerability$Month){
      for(j in Spatial_Vulnerability$Depth){
        
        Dens<-data.frame(Depth=DepthAll,Prob=scales::rescale(dlnorm(DepthAll,meanlog = Param[1],sdlog = Param[2],log = FALSE),to=c(0,1)))
        
        Spatial_Vulnerability$Vul[Spatial_Vulnerability$Month==i &Spatial_Vulnerability$Depth==j]<-Dens$Prob[which(Dens$Depth==j)]
      }
    }
  }
  
 return(Spatial_Vulnerability)
}

