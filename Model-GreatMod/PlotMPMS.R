library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(epimod)
library(ggh4x)
library(patchwork)

############################
##########  PLOT  ########## 
############################

load("MSMP1.Rdata")
load("MSMP2.Rdata")
load("MSMP2_1.Rdata")
load("MSMP2_2.Rdata")
load("MSMP3.Rdata")

trace = rbind(rbind(MSMP1$trace,MSMP2$trace),MSMP3$trace)
trace$ID = paste0("MSMP",trace$ID)
reference = rbind(rbind(MSMP1$ref,MSMP2$ref),MSMP3$ref)
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
  scale_color_manual(values = c("MSMP1"="#262035","MSMP2"= "#267578", "MSMP3"="#F1661D"))+
  theme(legend.position = "bottom",
        legend.box.background = element_rect(color="black", size=2),
        legend.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 12, angle = 45,hjust = 1)) +
  labs(x = "", y = "", col = "")

ggsave(pl2,filename = "TeffTreg.pdf",width = 10,height =8)
  
############################## 
### what if analysi on MSMP2

MSMP2_2$trace$ID = paste0(MSMP2_2$trace$ID," III injection")
MSMP2_1$trace$ID = "2"
trace2 = rbind(MSMP2_2$trace,MSMP2_1$trace)

vlines= data.frame(ATZv = c(0,360,1080),tyL=c("MPMS2","MPMS2","MPMS2 III injection"))
steps = c("BASELINE"," M6"," M12"," M18"," M24"," M36"," M48","M60"," M72" )
s = data.frame(step = c(0,360,1080),symbol = fa,Observation = 0,
               Place = rep(c("Antigen","Teff_out_t17","Teff_out_t1","Treg_out"),each = 3) )

reference$TID = paste0(reference$Time,"_",reference$ID)

cell_labels = c(`Teff_out_t1` = "Th1",`Teff_out_t17` = "Th17",`Treg_out` = "Treg")
pl = ggplot() 	+
  facet_wrap(~Place, scales = "free", labeller = as_labeller(cell_labels))+
  geom_text(data = s %>% filter(! Place %in% c("Teff_Treg","Antigen")),
            aes(x = step,y = Observation-1, label=symbol),nudge_x = 50,
            size = 4, family='fontawesome-webfont',col="#009E73")+
  geom_vline(data= vlines ,
             aes(xintercept = ATZv,linetype = tyL),
             col="#009E73",  alpha = .8, size = 1) +
  geom_line(data = trace2 %>% filter(! Place %in% c("Teff_Treg","Antigen")),
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

ggsave(pl,filename = "Plot/FigWIA.pdf",width = 14,height =10)

trace2$ID2 = paste0("W-IF ",trace2$ID)
trace2$ID4 = "What-if analysis MSMP2"
trace2$ID = gsub(trace2$ID,pattern = " III injection",replacement = "")
fract = rbind(rbind(MPMS1$trace, MPMS2$trace),MPMS3$trace)
fract$ID = paste0("MSMP",fract$ID)
fract$ID2 = "MSMPs"
fract$ID4 = "MSMPs"
fract = rbind(trace2,fract) %>% filter(Place == "Teff_Treg")
fract$ID3 = paste0(fract$ID,"_",fract$ID2)

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
  scale_linetype_manual(values = c("MSMPs"="solid","W-IF 2_2 III injection"= "dashed", "W-IF 2"="solid"))+
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

##########################
#### FIG3 of the paper ###
##########################

strip <- strip_themed(text_x = elem_list_text(colour  = c("#262035", "#267578", "#F1661D")))
pl_up = ggplot() 	+
  facet_grid2(~ID, scales = "free",strip = strip)+
  geom_text(data = s%>% filter(Place == "Teff_out_t17"), aes(x = step,y = Observation-1, label=symbol),nudge_x = 50,
            size = 4, family='fontawesome-webfont',col="#009E73")+
  geom_boxplot(data = reference %>% filter(Place == "Teff_out_t17"), aes(x = Time, y = Obs, group = Time),col="#5E5E5E",alpha = 0.5) +
  #geom_point(data = reference, aes(x = Time, y = Obs, group = IDref),col="#5E5E5E",alpha = 1)+
  geom_vline(data=vlines, aes(xintercept = ATZv ),
             col="#009E73", linetype = "solid", alpha = .8, size = 1) +
  geom_line(data = trace %>%filter(Place == "Teff_out_t17"),
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

pl_up_antig = ggplot() 	+
  facet_wrap(~Place, scales = "free")+
  geom_text(data = s%>% filter(Place == "Antigen"), aes(x = step,y = Observation-1, label=symbol),nudge_x = 50,
            size = 4, family='fontawesome-webfont',col="#009E73")+
  geom_vline(data=vlines, aes(xintercept = ATZv ),
             col="#009E73", linetype = "solid", alpha = .8, size = 1) +
  geom_line(data = trace %>%filter(Place == "Antigen", ID == "MSMP1"),
            aes(x = Time, y = Obs, group = ID),col = "black",linewidth = 1) +
  scale_x_continuous(breaks = unique(reference$Time),label = steps) +
  theme_bw()+
  theme(legend.position = "none",
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"),
        strip.text.y = element_text(
          size = 10, face = "bold.italic"),
        axis.title.y = element_text(face = "bold", size = 12),
        legend.box.background = element_rect(color="black", size=2),
        legend.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 12, angle = 45,hjust = 1)) +
  labs(x = "", y = "Cell quantification", col = "")

pl2 = ggplot() 	+
  facet_grid( ~ Place,labeller = as_labeller(c(`Teff_Treg` = "Teff / Treg")))+
  geom_vline(data=vlines, aes(xintercept = ATZv ),
             col="#009E73", linetype = "solid", alpha = .8, size = 1) +
  geom_line(data = trace%>%filter(Place == "Teff_Treg"),
            aes(x = Time, y = Obs,col = ID, group = ID),linewidth = 1) +
  scale_x_continuous(breaks = unique(reference$Time),label = steps) +
  theme_bw()+
  scale_color_manual(values = c("MSMP1"="#262035","MSMP2"= "#267578", "MSMP3"="#F1661D"))+
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
  labs(x = "", y = "", col = "")

pl3 = ggplot() 	+
  facet_wrap(~ID4,labeller = as_labeller(c(`What-if analysis MSMP2` = "Teff / Treg")) )+
  geom_text(data = s%>% filter(ID4 == "What-if analysis MSMP2"), aes(x = step,y = Observation-.1, label=symbol),nudge_x = 50,
            size = 4, family='fontawesome-webfont',col="#009E73")+
  geom_vline(data=vlines%>% filter(ID4 == "What-if analysis MSMP2"), aes(xintercept = ATZv),
             col="#009E73",  alpha = .8, size = 1) +
  geom_line(data = fract %>% filter(ID4 == "What-if analysis MSMP2") ,
            aes(x = Time, y = Obs, col = ID,  group = ID3, linetype = ID2),
            linewidth = 1,alpha=0.9) +
  scale_x_continuous(breaks = unique(reference$Time),label = steps) +
  scale_color_manual(values = c("2"="#267578","2_2"= "#267578"),label = c("W-IF 2_2 III injection", "W-IF 2"))+
  scale_linetype_manual(values = c("MSMPs"="solid","W-IF 2_2 III injection"= "dashed", "W-IF 2"="solid"))+
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
  labs(x = "", y = "Ratio - cell quantification", col = "",fill = "",linetype = "")


pl_bottom = ggplot() 	+
  facet_wrap(~Place, scales = "free", labeller = as_labeller(cell_labels))+
  geom_text(data = s %>% filter(! Place %in% c("Teff_Treg","Antigen")),
            aes(x = step,y = Observation-1, label=symbol),nudge_x = 50,
            size = 4, family='fontawesome-webfont',col="#009E73")+
  geom_vline(data= vlines ,
             aes(xintercept = ATZv,linetype = tyL),
             col="#009E73",  alpha = .8, size = 1) +
  geom_line(data = trace2 %>% filter(! Place %in% c("Teff_Treg","Antigen")),
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

pl_bottom_antig = ggplot() 	+
  facet_wrap(~Place, scales = "free")+
  geom_text(data = s %>% filter( Place %in% c("Antigen")),
            aes(x = step,y = Observation-1, label=symbol),nudge_x = 50,
            size = 4, family='fontawesome-webfont',col="#009E73")+
  geom_vline(data= vlines ,
             aes(xintercept = ATZv,linetype = tyL),
             col="#009E73",  alpha = .8, size = 1) +
  geom_line(data = trace2 %>% filter(Place %in% c("Antigen")),
            aes(x = Time, y = Obs, group = ID,linetype = ID),linewidth = 1,alpha=0.9) +
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


(pl_up_antig + labs(y="",title = "A.1)")) + (pl2+ labs(y="",title = "A.2)")) + (pl_up + labs(y="",title = "A.3)")) +
  (pl_bottom_antig + labs(y="",title = "B.1)")) + (pl3+ labs(y="",title = "B.2)")) + (pl_bottom + labs(y="",title = "B.3)")) +
  plot_layout(widths = c(1,1, 3)) & theme(legend.position = "none")

ggsave(,filename = "Plot/Fig3.pdf",width = 20,height =8)
