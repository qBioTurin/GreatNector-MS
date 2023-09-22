library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

plot.generation = function(folder,  reference = NULL){
  NumberMPMS=2
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
  
  trace <-traceDF
  reference <- referenceG
  trace$ID = paste0("MPMS",trace$ID)
  reference$IDref = reference$ID 
  reference$ID = paste0("MPMS",reference$ID)
  
  vlines= data.frame(ATZv = c(0,360),tyL="ATZ injections")
  steps = c("BASELINE"," M6"," M12"," M18"," M24"," M36"," M48","M60"," M72" )

  
  MPMS2_wi = list(reference,trace)
  
  pl = ggplot() 	+
    facet_wrap(~Place, scales = "free")+
    geom_boxplot(data = reference, aes(x = Time, y = Obs, group = Time),col="red",alpha = 0.5) +
    geom_point(data = reference, aes(x = Time, y = Obs, group = IDref),col="red",alpha = 1)+
    geom_vline(data=vlines, aes(xintercept = ATZv ),
               col="#009E73", linetype = "solid", alpha = .8, size = 1) +
    geom_line(data = trace%>%filter(Place != "Teff_Treg"),
              aes(x = Time, y = Obs, group = ID,col = ID),linewidth = 1) +
    scale_x_continuous(breaks = unique(reference$Time),label = steps) +
    scale_color_manual(values = c("MPMS1"="#262035","MPMS2"= "#267578", "MPMS3"="#F1661D"))+
    theme_bw()+
    theme(legend.position = "none",
          strip.text.x = element_text(
            size = 12, color = "white", face = "bold.italic"),
          strip.text.y = element_text(
            size = 10, face = "bold.italic"),
          legend.box.background = element_rect(color="black", size=2),
          legend.title = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold", size = 12),
          axis.text.x = element_text(face = "bold", size = 12, angle = 45,hjust = 1)) +
    labs(x = "", y = "", col = "")

return(list(pl = pl,trace = MPMS2_wi))
}
# 
# pl = plot.generation(folder = "MS_Model_calibration",
# 								reference = "./Input/referenceGroup2.csv",
# 								event_ATZ = c(180,540),
# 								event_Antig = as.integer(seq(30,2000,length.out = 4) ) )
# 
# 
# ggsave(pl$pl,filename= "calibration.pdf",device = "pdf",width = 12,height = 10)
# 




