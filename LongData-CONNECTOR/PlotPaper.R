library(dplyr)
library(readxl)
library(tidyr)
library(ggplot2)
library(patchwork)
TemporalCurve = readRDS(file = "Data/temporalCurveClustered.Rds")
MeanCurves = readRDS(file = "RData/DataWithoutCurves/MeanCurves.Rds")

steps = c("BASELINE"," M6"," M12"," M18"," M24"," M36"," M48","M60"," M72" )
names(steps) = c(0,6 ,12, 18, 24, 36, 48, 60, 72)

TemporalCurve$Step = steps[paste(TemporalCurve$Time)]
TemporalCurve$Step = factor(x = TemporalCurve$Step ,levels = steps)
TemporalCurve$Exp = gsub(TemporalCurve$Exp,pattern = "_prod",replacement = "")
TemporalCurve = TemporalCurve %>% 
  mutate( Discrete.AGE.AT.ONSET =  ifelse(AGE.AT.ONSET<30, "under 30", "over equal 30"),
          Discrete.DISEASE.DURATION =  ifelse(DISEASE.DURATION<5, "less 5", "greater equal 5") )

MeanCurves$Step = steps[paste(MeanCurves$Time)]
MeanCurves$Step = factor(x = MeanCurves$Step ,levels = steps)
MeanCurves$Exp = gsub(MeanCurves$Exp,pattern = "_prod",replacement = "")

TemporalCurve$Exp = factor(x = TemporalCurve$Exp ,levels = c("CD4","Th1","Th17","Treg","EDSS"))
MeanCurves$Exp = factor(x = MeanCurves$Exp ,levels = c("CD4","Th1","Th17","Treg","EDSS"))

EDSSCurve = TemporalCurve %>%
  filter(Exp == "EDSS") %>%
  dplyr::select(ID, Cluster) %>% 
  distinct() %>%
  rename(EDSS = Cluster)

TemporalCurve = TemporalCurve %>% filter(Exp != "EDSS")
TemporalCurve = merge(TemporalCurve,EDSSCurve,by="ID",all = T)
MeanCurves = MeanCurves %>% filter(Exp != "EDSS")


rel = grep(pattern = "RELAPSES", colnames(TemporalCurve))
do.call(rbind,
        lapply(TemporalCurve %>% 
                 filter(Relapses!="No rel") %>% 
                 select(ID) %>%
                 distinct() %>% 
                 unlist(),
               function(i){
                 TemporalCurve_sub = TemporalCurve %>% filter(ID == i)
                 ftime = gsub(pattern = "M",
                              replacement = "",
                              x = TemporalCurve_sub[1,rel][!is.na(TemporalCurve_sub[1,rel] )]
                 )
                 ftime = unlist( strsplit(x = ftime,split = ",") )
                 ftime = min(as.numeric(ftime))
                 TemporalCurve_sub = TemporalCurve_sub %>% 
                   dplyr::select(ID,Time, Exp) %>% 
                   na.omit() %>%
                   dplyr::group_by(Exp) %>%
                   dplyr::mutate(dif = Time - ftime ) %>% 
                   filter(dif > 0)  %>% 
                   filter(dif == min(dif)) %>% 
                   ungroup() %>%
                   dplyr::select(-dif)
                 TemporalCurve_sub$FirstRel = ftime
                 return(TemporalCurve_sub)
               })
) -> firstRelapse
firstRelapse = merge(firstRelapse,TemporalCurve,all.x = T,by = c("ID","Time", "Exp")) 
firstRelapse$Exp[firstRelapse$ID =="GA0001"] = "Treg"

featCol = c("ID","GENDER","Discrete.AGE.AT.ONSET","Discrete.DISEASE.DURATION","PREVIOUS.IMMUNOMODULANT.THERAPY","PREVIOUS.IMMUNOSUPPRESSANT.THERAPY","Aut..post.ATZ","Relapses","EDSS")

tmp = TemporalCurve%>%
  arrange(Time) %>% 
  drop_na(Observation)

tmp = tmp %>% mutate(Relapses = ifelse(Relapses=="No rel","Zero","At least one"))
tmp = tmp %>% mutate(ClusterKmeans = paste0("MSMP",ClusterKmeans))
firstRelapse = firstRelapse %>% mutate(ClusterKmeans = paste0("MSMP",ClusterKmeans))

load.fontawesome()
search_emoji('syring')
fa <- fontawesome(c('fa-bolt'))
s1 = tmp %>% group_by(Exp) %>% filter(Observation == max(Observation)) %>% select(Observation,Exp)
s1$symbol = fa
s1$step = "BASELINE"
s2 = tmp %>% group_by(Exp) %>% filter(Observation == max(Observation)) %>% select(Observation,Exp)
s2$symbol = fa
s2$step = " M12"
s = rbind(s1,s2)

vlines= data.frame(ATZv = c("BASELINE"," M12"),tyL="ATZ injections")

### Small figures for Fig.1

for( i in unique(tmp$Exp)){
  pltmp = tmp %>% filter(Exp == i) %>%
    ggplot() +
    geom_line(aes(x = Step, y = Observation,group = ID, col = ID),alpha = 0.6, size = 1.5)+
    geom_line(data = MeanCurves%>% filter(Exp == i), aes(x = Step, y = Mean, group = Exp),
              linetype = "dashed", size =1)+
    facet_wrap( ~Cluster,scales = "free",ncol = 1) +
    theme_bw()+
    scale_x_discrete(breaks = steps,limit = steps)+
    theme(legend.position = "none",
          legend.box.background = element_rect(color="black", size=2),
          legend.title = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold", size = 12),
          axis.text.x = element_text(face = "bold", size = 12, angle = 45,hjust = 1)) +
    labs(x="",y="",col = "")
  assign(paste0("pl",i),pltmp)
}
plCD4 + plTh1 + plTh17 + plTreg

ggsave(,filename = "Fig1curves.pdf",device = "pdf",width = 16,height = 16)

### FIGURA 2A
ggplot(tmp) +
  geom_line(aes(x = Step, y = Observation,group = ID, col = Relapses),size = 1.5, alpha = 0.6)+
  facet_wrap( ~Exp,scales = "free_y",ncol = 1) +
  theme_bw()+
  scale_x_discrete(breaks = steps,limit = steps)+
  geom_text(data = s, aes(x = step,y = 0, label=symbol),
            nudge_x = -0.08,
            size = 4, family='fontawesome-webfont',col="#009E73")+
  theme(legend.position = "bottom",
        legend.box.background = element_rect(color="black", size=2),
        legend.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 12, angle = 45,hjust = 1)) +
  geom_vline(data=vlines, aes(xintercept = ATZv ),col="#009E73", linetype = "solid", alpha = .8, size = 1) +
  labs(x="",y="",col = "Number of\nrelapses")+
  scale_color_manual(values = c("#FC4E07", "#56B4E9"))

ggsave(,filename = "Fig2A.pdf",device = "pdf",width = 6,height = 10)

meanC = tmp %>% 
  group_by(Exp,ClusterKmeans,Step) %>%
  summarize(M = mean(Observation,na.rm = T)) %>%
  ungroup()%>%
  as.data.frame()


ggplot() +
  geom_text(data = s, aes(x = step,y = 0, label=symbol),
            nudge_x = -0.08,
            size = 4, family='fontawesome-webfont',col="#009E73")+
  geom_line(data = tmp %>% filter(Relapses == "Zero"),
            aes(x = Step, y = Observation,group = ID), col = "#56B4E9",alpha = 0.7,size =1)+
  geom_line(data = tmp %>% filter(Relapses != "Zero") %>% mutate(ID = paste0(ClusterKmeans,"-",ID)),
            aes(x = Step, y = Observation,group = ID, col = ID),alpha = 0.8,size =1.5)+
 geom_point(data = firstRelapse %>%
              filter(Relapses != "Zero", Exp == "Treg") %>%
              mutate(ID = paste0(ClusterKmeans,"-",ID),
                     Obs = -5) %>% 
              dplyr::select(Cluster,TypeRel,ID, Time,Obs,  Exp, FirstRel,ClusterKmeans),
            aes(x = 1+FirstRel/9, y = Obs, fill = ID,col = ID),
            size = 5,shape = 25 )+
  geom_line(data = meanC,aes(x = Step, y = M,group =Exp),color="black",linetype = "dashed",size = 1)+
  facet_grid( Exp~ClusterKmeans,scales = "free_y") +
  theme_bw()+
  scale_x_discrete(breaks = steps,limit = steps)+
  expand_limits(y=0)+
  coord_cartesian( clip="off")+
  guides(col=guide_legend(nrow=3,ncol = 3,byrow=TRUE),
         fill=guide_legend(nrow=3,ncol = 3,byrow=TRUE))+
  theme(legend.position = "bottom",plot.margin = unit(c(1,1,1,0), "lines"),
        legend.box.background = element_rect(color="black", size=2),
        legend.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 12, angle = 45,hjust = 1)) +
  geom_vline(data=vlines, aes(xintercept = ATZv ),col="#009E73", linetype = "solid", alpha = .7, size = 1) +
  labs(x="",y="",col = "Number of\nrelapses",fill ="Number of\nrelapses" )+
  scale_color_manual(values = rep(c("#FC4E07", "#E69F00", "#CC79A7"),3))+
scale_fill_manual(values = rep(c("#FC4E07", "#E69F00", "#CC79A7"),3))

ggsave(,filename = "Fig2B_mean.pdf",device = "pdf",width = 14,height = 10)


ggplot(tmp) +
  geom_line(aes(x = Step, y = Observation,group = ID, col = Relapses),alpha = 0.6)+
  geom_point(data = firstRelapse,
             aes(x = Step, y = Observation, group = ID, col = Relapses),size = 2 )+
  facet_grid( Exp~ClusterKmeans,scales = "free_y") +
  theme_bw()+
  scale_x_discrete(breaks = steps,limit = steps)+
  theme(legend.position = "none",
        legend.box.background = element_rect(color="black", size=2),
        legend.title = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold", size = 12, angle = 45,hjust = 1)) +
  geom_vline(data=vlines, aes(xintercept = ATZv ),col="#009E73", linetype = "dashed", alpha = .9, size = 1) +
  labs(x="",y="",col = "Number of\nrelapses")+
  scale_color_manual(values = c("#FC4E07", "#56B4E9"))

ggsave(,filename = "Fig2B.pdf",device = "pdf",width = 11,height = 10)
