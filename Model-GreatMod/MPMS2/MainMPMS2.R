# library(devtools)
# install_github("https://github.com/qBioTurin/epimod")
# downloadContainers(tag = "latest")

library(epimod)
setwd("~/GreatNector-MS/Model-GreatMod/MPMS2/")
source("../RFunctions/PlotFig4.R")

run = function(init,solver = "MS_Model.solver",ATZinj = c(180, 540)){

    model.analysis(solver_fname =solver ,
                   f_time = 84*30, # days in 84 months
                   s_time = 1,
                   parameters_fname = "../Input/ParamsList2.csv",
                   functions_fname = "../RFunctions/Functions.R",
                   ini_v = init, # debug = T,
                   event_function = "EventFun",
                   event_times = ATZinj
    )

  pl = plot.generation(folder = paste0(gsub(pattern = ".solver",replacement = "",x = solver),"_analysis"),
                       reference = "../Input/referenceGroup2.csv" )
  
  
  return(pl)
}

# 1) Generation of the model starting from its graphical representation

model.generation(net_fname = "../Net/MS_Model.PNPRO",
								 transitions_fname = "../Net/GenTransitions.cpp")

## If you want to run the what-if analysis then the solver to exploit is MS_Model2.solver

model.generation(net_fname = "../Net/MS_Model.PNPRO",
                 transitions_fname = "../Net/GenTransitions_2.cpp")
file.rename("MS_Model.solver","MS_Model2.solver")


# 2) Model Analysis

parmNamesMS = c("TeDup","TrDup2","TrDup","TeDup2",
                "TrkTe","TekA","Pass_BBB","ATZkill",
                 "VentryTreg","VentryTeff_1", "VentryTeff_17",
                "Consuption","A","B","Period")

saveRDS(parmNamesMS,file = "../Input/params_name2.RDS")
init = rep(0,length(parmNamesMS))
names(init) = parmNamesMS

init["TeDup"] = .9998
init["TeDup2"] = 0.00024 # dip from antigen
init["TrDup2"] = 0.00015  # dip from teff
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

run1 = run(init,solver = "MS_Model2.solver")

##### WHAT IF analysis plot:

# III injection at 36 month
run2 = run(init,solver = "MS_Model2.solver",ATZinj = 180+c(0,360,1080)) 

##### plot
run2$trace[[2]]$ID = paste0(run2$trace[[2]]$ID," III injection")
run2$trace[[1]]$ID = paste0(run2$trace[[1]]$ID," III injection")
trace = rbind(run1$trace[[2]],run2$trace[[2]])

vlines= data.frame(ATZv = c(0,360,1080),tyL=c("MPMS2","MPMS2","MPMS2 III injection"))
steps = c("BASELINE"," M6"," M12"," M18"," M24"," M36"," M48","M60"," M72" )

library(emojifont)
load.fontawesome()
fa <- fontawesome(c('fa-bolt'))
s = data.frame(step = c(0,360,1080),symbol = fa,Observation = 0,
               Place = rep(c("Antigen","Teff_out_t17","Teff_out_t1","Treg_out"),each = 3) )
load("../MPMS1.Rdata")
load("../MPMS2.Rdata")
load("../MPMS3.Rdata")

reference = rbind(rbind(MPMS1$ref, MPMS2$ref),MPMS3$ref)
reference$ID = as.factor(paste0("MSMP",reference$ID))
reference$TID = paste0(reference$Time,"_",reference$ID)

pl = ggplot() 	+
  facet_wrap(~Place, scales = "free")+
  geom_text(data = s, aes(x = step,y = Observation-1, label=symbol),nudge_x = 50,
            size = 4, family='fontawesome-webfont',col="#009E73")+
  geom_vline(data=vlines, aes(xintercept = ATZv,linetype = tyL),
             col="#009E73",  alpha = .8, size = 1) +
  geom_line(data = trace%>%filter(Place != "Teff_Treg"),
            aes(x = Time, y = Obs, group = ID,linetype = ID),linewidth = 1,alpha=0.9) +
  geom_boxplot(data = reference, aes(x = Time, y = Obs, group = TID,fill = ID, col = ID),alpha=0.6) +
  scale_x_continuous(breaks = unique(reference$Time),label = steps) +
  scale_color_manual(values = c("MSMP1"="#262035","MSMP2"= "#267578", "MSMP3"="#F1661D"))+
  scale_fill_manual(values = c("MSMP1"="#262035","MSMP2"= "#267578", "MSMP3"="#F1661D"))+
  theme_bw()+
  guides(linetype = "none")+
  theme(legend.position = "bottom",
        #strip.background.x = element_rect(fill ="#267578"),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"
          ),
        axis.title.y = element_text(face = "bold", size = 12),
        legend.box.background = element_rect(color="black", size=2),
        legend.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 12, angle = 45,hjust = 1)) +
  labs(x = "", y = "Cell quantification", col = "",fill = "")
pl
ggsave(pl,filename = "FigWIA.pdf",width = 14,height =10)

trace$ID2 = paste0("W-IF ",trace$ID)
trace$ID4 = "What-if analysis MSMP2"
trace$ID = gsub(trace$ID,pattern = " III injection",replacement = "")
fract = rbind(rbind(MPMS1$trace, MPMS2$trace),MPMS3$trace)
fract$ID = paste0("MSMP",fract$ID)
fract$ID2 = "MSMPs"
fract$ID4 = "MSMPs"
fract = rbind(trace,fract) %>% filter(Place == "Teff_Treg")
fract$ID3 = paste0(fract$ID,"_",fract$ID2)
fract$ID = gsub(fract$ID,pattern = "MPMS",replacement = "MSMP")
fract$ID2 = gsub(fract$ID2,pattern = "MPMS",replacement = "MSMP")
fract$ID3 = gsub(fract$ID3,pattern = "MPMS",replacement = "MSMP")
s = rbind(s%>%filter(step!=1080)%>% mutate(ID4 = "MSMPs"),s%>% mutate(ID4 = "What-if analysis MSMP2"))
vlines = rbind(vlines%>%filter(ATZv!=1080)%>% mutate(ID4 = "MSMPs"),vlines%>% mutate(ID4 = "What-if analysis MSMP2"))
pl2 = ggplot() 	+
  facet_wrap(~ID4,nrow = 1 )+
  geom_text(data = s, aes(x = step,y = Observation-.1, label=symbol),nudge_x = 50,
            size = 4, family='fontawesome-webfont',col="#009E73")+
  geom_vline(data=vlines, aes(xintercept = ATZv),
             col="#009E73",  alpha = .8, size = 1) +
  geom_line(data = fract,
            aes(x = Time, y = Obs, col = ID,  group = ID3, linetype = ID2),
            linewidth = 1,alpha=0.9) +
  scale_x_continuous(breaks = unique(reference$Time),label = steps) +
  scale_color_manual(values = c("MSMP1"="#262035","MSMP2"= "#267578", "MSMP3"="#F1661D"))+
  scale_linetype_manual(values = c("MSMPs"="solid","W-IF MSMP2 III injection"= "dashed", "W-IF MSMP2"="solid"))+
  theme_bw()+
  theme(legend.position = "bottom",
        #strip.background.x = element_rect(fill ="#267578"),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"
        ),
        legend.box.background = element_rect(color="black", size=2),
        legend.title = element_text(face = "bold"),
       axis.title.y = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 12, angle = 45,hjust = 1)) +
  labs(x = "", y = "Teff / Treg", col = "",fill = "",linetype = "")
ggsave(pl2,filename = "FigRapp.pdf",width = 10,height =6)




