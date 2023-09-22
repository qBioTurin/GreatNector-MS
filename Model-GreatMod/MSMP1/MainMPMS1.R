# library(devtools)
# install_github("https://github.com/qBioTurin/epimod", ref="variabilityFBA")
# downloadContainers(tag = "2.0.0")

library(epimod)
setwd("MPMS1/")
source("./RFunctions/Plot.R")

run = function(init){

    model.analysis(solver_fname = "Net/MS_Model.solver",
                   f_time = 84*30, # days in 84 months
                   s_time = 1,
                   parameters_fname = "Input/ParamsList.csv",
                   functions_fname = "RFunctions/Functions.R",
                   ini_v = init, # debug = T,
                   event_function = "EventFun",
                   event_times = sort(c(180,540 ))
    )
  
  
  pl = plot.generation(folder = "MS_Model_analysis",
                       reference = "Input/referenceGroup1.csv",
                       event_ATZ = c(0,360), # Idose = after 6M and the IIdose after 12M w.r.t. the Idose
                       event_Antig = 0 )
  return(pl$pl)
}

# 1) Generation of the model starting from its graphical representation

model.generation(net_fname = "Net/MS_Model.PNPRO",
								 transitions_fname = "Net/GenTransitions.cpp")
system("mv MS_Model.* ./Net")

# 2) Model Analysis

parmNamesMS = c("TeDup","TrDup2","TrDup","TeDup2",
                "TrkTe","TekA","Pass_BBB","ATZkill",
                 "VentryTreg","VentryTeff_1", "VentryTeff_17",
                "Consuption","A","B","Period")

saveRDS(parmNamesMS,file = "./Input/params_name.RDS")
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

run(init)+labs(title = "MPMS1")

