library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(epimod)
library(ggh4x)

trace.generation = function(folder,NumberMPMS, reference){
  
  file <- list.files(paste0("./",folder,"/"),pattern='*.trace')

  trace <- data.frame(ID = NumberMPMS, read.csv(paste0('./',folder,'/',file), sep="",stringsAsFactors = F))
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
    referenceG$ID = NumberMPMS

  return(list(ref = referenceG,trace = traceDF))
}
run = function(init,NumberMPMS){
  setwd(paste0("./MPMS",NumberMPMS))
  
  model.analysis(solver_fname = "Net/MS_Model.solver",
                 f_time = 84*30, # days in 84 months
                 s_time = 1,
                 parameters_fname = "Input/ParamsList.csv",
                 functions_fname = "RFunctions/Functions.R",
                 ini_v = init, # debug = T,
                 event_function = "EventFun",
                 event_times = sort(c(180,540 ))
  )
  
  
  tr = trace.generation(folder = "MS_Model_analysis",
                       reference = paste0("Input/referenceGroup",NumberMPMS,".csv"),
                       NumberMPMS = NumberMPMS)
  setwd("../")
  
  return(tr)
}
parmNamesMS = c("TeDup","TrDup2","TrDup","TeDup2",
                "TrkTe","TekA","Pass_BBB","ATZkill",
                "VentryTreg","VentryTeff_1", "VentryTeff_17",
                "Consuption")
init = rep(0,length(parmNamesMS))
names(init) = parmNamesMS


##########  MSMP1 ########## 
############################
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

MPMS1 = run(init,NumberMPMS = 1)

############################
##########  MSMP2 ########## 
############################

MPMS2 = run(init,NumberMPMS = 2)

############################
##########  MSMP3 ########## 
############################

MPMS3 = run(init,NumberMPMS = 3)

############################
##########  PLOT  ########## 
############################

save(MPMS1,file = "MPMS1.Rdata")
save(MPMS2,file = "MPMS2.Rdata")
save(MPMS3,file = "MPMS3.Rdata")
trace = rbind(rbind(MPMS1$trace,MPMS2$trace),MPMS3$trace)
trace$ID = paste0("MSMP",trace$ID)
reference = rbind(rbind(MPMS1$ref,MPMS2$ref),MPMS3$ref)
reference$ID = paste0("MSMP",reference$ID)

vlines= data.frame(ATZv = c(0,360),tyL="ATZ injections")
steps = c("BASELINE"," M6"," M12"," M18"," M24"," M36"," M48","M60"," M72" )
strip <- strip_themed(background_x = elem_list_rect(fill = c("#262035", "#267578", "#F1661D")))
library(emojifont)
load.fontawesome()
fa <- fontawesome(c('fa-bolt'))
s = data.frame(step = c(0,360),symbol = fa,Observation = 0,
               Place = rep(c("Antigen","Teff_out_t17","Teff_out_t1","Treg_out"),each = 2) )

pl = ggplot() 	+
  facet_grid2(Place~ID, scales = "free",strip = strip)+
  geom_text(data = s, aes(x = step,y = Observation-1, label=symbol),nudge_x = 50,
            size = 4, family='fontawesome-webfont',col="#009E73")+
  geom_boxplot(data = reference, aes(x = Time, y = Obs, group = Time),col="#5E5E5E",alpha = 0.5) +
  #geom_point(data = reference, aes(x = Time, y = Obs, group = IDref),col="#5E5E5E",alpha = 1)+
  geom_vline(data=vlines, aes(xintercept = ATZv ),
           col="#009E73", linetype = "solid", alpha = .8, size = 1) +
  geom_line(data = trace%>%filter(Place != "Teff_Treg"),
            aes(x = Time, y = Obs, group = ID,col = ID),linewidth = 1) +
  scale_x_continuous(breaks = unique(reference$Time),label = steps) +
  scale_color_manual(values = c("MSMP1"="#262035","MSMP2"= "#267578", "MSMP3"="#F1661D"))+
  theme_bw()+
  theme(legend.position = "none",
        strip.text.x = element_text(
          size = 12, color = "white", face = "bold.italic"),
        strip.text.y = element_text(
          size = 10, face = "bold.italic"),
        axis.title.y = element_text(face = "bold", size = 12),
        legend.box.background = element_rect(color="black", size=2),
        legend.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 12, angle = 45,hjust = 1)) +
  labs(x = "", y = "Cell quantification", col = "")

pl
ggsave(pl,filename = "Fig3.pdf",width = 14,height =10)

pl2 = ggplot() 	+
  geom_vline(data=vlines, aes(xintercept = ATZv ),
             col="#009E73", linetype = "solid", alpha = .8, size = 1) +
  geom_line(data = trace%>%filter(Place == "Teff_Treg"),
            aes(x = Time, y = Obs,col = ID, group = ID),linewidth = 1) +
  scale_x_continuous(breaks = unique(reference$Time),label = steps) +
  theme_bw()+
  scale_color_manual(values = c("MPMS1"="#262035","MPMS2"= "#267578", "MPMS3"="#F1661D"))+
  theme(legend.position = "bottom",
        legend.box.background = element_rect(color="black", size=2),
        legend.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 12, angle = 45,hjust = 1)) +
  labs(x = "", y = "", col = "")

ggsave(pl2,filename = "TeffTreg.pdf",width = 10,height =8)
  

