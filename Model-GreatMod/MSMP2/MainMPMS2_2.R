# library(devtools)
# install_github("https://github.com/qBioTurin/epimod", ref="variabilityFBA")
# downloadContainers(tag = "2.0.0")

library(epimod)
setwd("./temporalanalysis_MS/MPMS2/")
source("./RFunctions/PlotFig4.R")

run = function(init,solver = "MS_Model.solver",ATZinj = c(180, 540)){

    model.analysis(solver_fname =paste0("Net/",solver),
                   f_time = 84*30, # days in 84 months
                   s_time = 1,
                   parameters_fname = "Input/ParamsList2.csv",
                   functions_fname = "RFunctions/Functions.R",
                   ini_v = init, # debug = T,
                   event_function = "EventFun",
                   event_times = ATZinj
    )

  pl = plot.generation(folder = paste0(gsub(pattern = ".solver",replacement = "",x = solver),"_analysis"),
                       reference = "Input/referenceGroup2.csv" )
  
  
  return(pl)
}

# 1) Generation of the model starting from its graphical representation

# first model
model.generation(net_fname = "Net/MS_Model.PNPRO",
                 transitions_fname = "Net/GenTransitions.cpp")
system("mv MS_Model.* ./Net")

## second model
model.generation(net_fname = "Net/MS_Model.PNPRO",
								 transitions_fname = "Net/GenTransitions_2.cpp")
file.rename("MS_Model.solver","MS_Model2.solver")
file.rename("MS_Model.net","MS_Model2.net")
file.rename("MS_Model.def","MS_Model2.def")
system("mv MS_Model2.* ./Net")
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
init["A"] = 5
init["B"] = 4.9
init["Period"] = 41

run1 = run(init,solver = "MS_Model.solver")

# III injection at 36 month
run2 = run(init,solver = "MS_Model2.solver",ATZinj = 180+c(0,360,1080)) 
run2$trace[[2]]$ID = paste0(run2$trace[[2]]$ID," III injection")
run2$trace[[1]]$ID = paste0(run2$trace[[1]]$ID," III injection")

