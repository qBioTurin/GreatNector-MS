library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(epimod)
library(ggh4x)
library(patchwork)

trace.generation = function(folder,NumberMSMP, reference){
  
  file <- list.files(paste0("./",folder,"/"),pattern='*.trace')
  
  trace <- data.frame(ID = NumberMSMP, read.csv(paste0('./',folder,'/',file), sep="",stringsAsFactors = F))
  traceDF =	trace %>% 
    mutate(Teff_Treg = Teff_in/(Treg_in + 1e-6) )%>%
    dplyr::select(-Teff_in,-Treg_in,-ATZ) %>%
    gather(Place, Obs,-ID, -Time) 
  
  traceDF = traceDF %>% mutate(Time = Time - 30*6)
  
  
  steps = seq(min(traceDF$Time),max(c(traceDF$Time)),by = 30*6)
  
  ### reference
  referenceG <- as.data.frame(read.csv(reference,
                                       header = FALSE,
                                       sep = "",stringsAsFactors = F))
  colnames(referenceG) = c("Time","ID","Place","Obs")
  
  referenceG[which(referenceG$Place == "Th1"),"Place"]  = "Teff_out_t1"
  referenceG[which(referenceG$Place == "Th17"),"Place"]  = "Teff_out_t17"
  referenceG[which(referenceG$Place == "Treg"),"Place"]  = "Treg_out"
  
  referenceG = referenceG %>%
    filter(Time <= max(traceDF$Time)) %>% 
    filter(Place != "ATZ")
  referenceG$IDref = referenceG$ID
  referenceG$ID = NumberMSMP
  
  return(list(ref = referenceG,trace = traceDF))
}
run = function(init,NumberMSMP){
  if(NumberMSMP == "2_2"){
    setwd("./MSMP2")
    solver = "Net/MS_Model2.solver"
    ATZinj = 180+c(0,360,1080)
    NumberMSMP2 = "2"
  }else if(NumberMSMP == "2_1"){
    setwd("./MSMP2")
    solver = "Net/MS_Model2.solver"
    ATZinj = 180+c(0,360)
    NumberMSMP2 = "2"
  }else{
    setwd(paste0("./MSMP",NumberMSMP))
    solver = "Net/MS_Model.solver"
    ATZinj = c(180,540 )
    NumberMSMP2 = NumberMSMP
  }
  
  
  model.analysis(solver_fname = solver,
                 f_time = 84*30, # days in 84 months
                 s_time = 1,
                 parameters_fname = "Input/ParamsList.csv",
                 functions_fname = "RFunctions/Functions.R",
                 ini_v = init, # debug = T,
                 event_function = "EventFun",
                 event_times = ATZinj
  )
  
  tr = trace.generation(folder = paste0(gsub(pattern = ".solver|Net",replacement = "",
                                             x = solver),"_analysis"),
                        reference = paste0("Input/referenceGroup",NumberMSMP2,".csv"),
                        NumberMSMP = NumberMSMP)
  setwd("../")
  
  return(tr)
}


parmNamesMS = c("TeDup","TrDup2","TrDup","TeDup2",
                "TrkTe","TekA","Pass_BBB","ATZkill",
                "VentryTreg","VentryTeff_1", "VentryTeff_17",
                "Consuption","A","B","Period")

init = rep(0,length(parmNamesMS))
names(init) = parmNamesMS

init["TeDup"] = .9998
init["TeDup2"] = 0.00024 # quest dip da antigen
init["TrDup2"] = 0.00015  # quest dip da teff
init["TrDup"] = 0.47
init["TrkTe"] = 0.01
init["TekA"] = 0.09      
init["Pass_BBB"] = 1e-4
init["ATZkill"] = 4e-02
init["VentryTreg"] = 0.72
init["VentryTeff_1"] = 0.82
init["VentryTeff_17"] = 0.11
init["Consuption"] = 0.06
init["A"] = 1
init["B"] = 1
init["Period"] = 180

##########  MSMP1 ########## 
############################


MSMP1 = run(init,NumberMSMP = 1)

############################
##########  MSMP2 ########## 
############################

MSMP2 = run(init,NumberMSMP = 2)

############################
##########  MSMP3 ########## 
############################

MSMP3 = run(init,NumberMSMP = 3)

save(MSMP1,file = "MSMP1.Rdata")
save(MSMP2,file = "MSMP2.Rdata")
save(MSMP3,file = "MSMP3.Rdata")

############################################ 
##########  MSMP2 WhatIf analysis ########## 
############################################ 
init2 = init
init2["A"] = 5
init2["B"] = 4.9
init2["Period"] = 41

MSMP2_1 = run(init2,NumberMSMP = "2_1")
MSMP2_2 = run(init2,NumberMSMP = "2_2")

save(MSMP2_1,file = "MSMP2_1.Rdata")
save(MSMP2_2,file = "MSMP2_2.Rdata")
