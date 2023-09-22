library(readr)
library(ggplot2)
library(ggthemes)
setwd("Results/")
SubEntities = c("Treg_out","IL10_out", "BBB",
                "Teff_in","IL17_in","IFNg_in","ODC_le1")

Entities = c("Teff_out","Treg_out","Antigen" ,
             "EffectorMemory","NK_out","BBB",
             "IL10_out","IL17_out","IFNg_out", 
             "Treg_in", "Teff_in","ODC_le1" ,
             "IL17_in","IFNg_in","IL10_in")

DataExcel <- read_table2("../input/DataExcel.csv")

ref <-lapply(2:(length(DataExcel[1,])), function(i){
  if(names(DataExcel[,i]) %in% c("IFNg_out", "IL17_out" ,"IL10_out", "Teff_out","Treg_out") ) t = 2+18/24
  else t = 4+18/24
  a<- data.frame(DataExcel[,c(1,i)], Entity = names(DataExcel[,i]), time = t )
  colnames(a) <- c("ID", "ref" ,"Entity","Time")
  return(a)
})

ref <- do.call("rbind",ref)

sub_entity = unique(ref$Entity)


### reading deterministic solution

traceH_Det <-data.frame(read.csv('./Healthy_model_analysis/MS_Model-analysys-1.trace', sep=""), patient = "Healthy")
traceS_Det <-data.frame(read.csv('./MS_model_analysis/MS_Model-analysys-1.trace', sep=""), patient = "MS")
trace_Det <-rbind(traceH_Det,traceS_Det)

DatiPlotlist_Det<- lapply(1:length(Entities),function(i){
  place<-Entities[i]
  
  sub_dati=trace_Det[,c("Time",place)]
  
  if(place == "ODC_le1") place = "ODC irr. damaged"
  dt=data.frame(sub_dati,ID=trace_Det$patient,Ent=place)
  colnames(dt)=c("Time","Value","ID","Entity")
  return(dt)
})
DatiPlot_Det<-do.call("rbind", DatiPlotlist_Det)

### reading stochastic simulations:
traceH_Stoch <-data.frame(read.csv('./Healthy_Stoch_analysis/MS_Model-analysys-1.trace', sep=""), patient = "Healthy")
traceS_Stoch <-data.frame(read.csv('./MS_Stoch_analysis/MS_Model-analysys-1.trace', sep=""), patient = "MS")

trace_Stoch <- rbind(traceS_Stoch,traceH_Stoch)

time <- unique(trace_Stoch$Time)
nsim <- unique(table(trace_Stoch$Time))

if (length(nsim)>1) cat("error")

trace_Stoch$Id = rep(1:nsim,each = length(time))

MeanData <-aggregate(trace_Stoch[, -which(colnames(trace_Stoch) %in% c("Time","Id","patient"))], by = list(trace_Stoch$Time,trace_Stoch$patient), mean)

DatiPlotlist_Stoch<- lapply(1:length(Entities),function(i){
  place<-Entities[i]
  
  
  sub_dati=trace_Stoch[,c("Time",place)]
  if(place == "ODC_le1") place = "ODC irr. damaged"
  dt=data.frame(sub_dati,ID=trace_Stoch$patient,SimID=trace_Stoch$Id,Ent=place)
  colnames(dt)=c("Time","Value","ID","SimID","Entity")
  return(dt)
})

MeanDatalist<- lapply(1:length(Entities),function(i){
  place<-Entities[i]
  
  sub_dati=MeanData[,c(colnames(MeanData)[1],place)]
  if(place == "ODC_le1") place = "ODC irr. damaged"
  dt=data.frame(sub_dati,ID=MeanData[,2],Ent=place)
  
  colnames(dt)=c("Time","Value","ID","Entity")
  return(dt)
})

DatiPlot<-do.call("rbind", DatiPlotlist_Stoch)
Mean <- do.call("rbind", MeanDatalist)


Plotgeneration<- function(DataFrame, Stoch, Mean=NULL, Entities)
{
  DatiPlot = DataFrame[which(DataFrame$Entity %in% Entities),]
  DatiPlot$Entity <- factor(DatiPlot$Entity, levels = Entities)
  
  ref = ref[which(ref$Entity %in% Entities),]
  
  pl<- ggplot() 
  
  if(Stoch){
    MeanDataF <- Mean[which(Mean$Entity %in% Entities ),]
    
    pl <- pl + geom_line(data = DatiPlot, aes(x=Time/24,y=Value,group = SimID,  col="Simulations"),alpha=.4)+
      geom_line(data = MeanDataF ,aes(x=Time/24,y=Value,col = ID))+
      scale_colour_manual("",values = c("Simulations" = "grey" ,"MS"="red","Healthy"="blue"), labels = c("Healthy","MS","Simulations") )+
      facet_grid(Entity ~ ID, scales = "free") + 
      geom_vline(data = ref, aes( col= ID , xintercept = Time), linetype = "longdash" )
  }else{
    pl <- pl + geom_line(data = DatiPlot, aes(x=Time/24,y=Value, group = ID, col= ID))+
      scale_colour_manual("",values = c("MS"="red","Healthy"="blue"), labels = c("Healthy","MS") )+
      facet_wrap(Entity~., scales = "free",ncol=3)+
      geom_vline(data = ref, aes( xintercept = Time), linetype = "longdash", col = "gray" )
  }
  
  pl = pl + geom_violin(data = ref, aes(x = Time, y = ref,fill=ID, col= ID ),alpha=.5)+
    theme(axis.text=element_text(size = 25, hjust = 0.5),
          axis.text.x=element_text(angle=+90,vjust=0.5, hjust=1),
          axis.title=element_text(size=28,face="bold"),
          axis.line = element_line(colour="black"),
          plot.title=element_text(size=30, face="bold", vjust=1, lineheight=0.6),
          legend.text=element_text(size=24),
          legend.position="bottom",
          legend.title=element_blank(),
          legend.key=element_blank(),
          legend.key.size = unit(.9, "cm"),
          legend.key.width = unit(.9,"cm"),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          plot.margin=unit(c(5,5,5,5),"mm"),
          strip.text = element_text(size = 20))+
    labs(x="Days", y="Quantity" )+
    scale_fill_manual("",values = c("MS"="red","Healthy"="blue","Simulations" = NULL),
                      labels = c("Healthy","MS","Simulations")  )
  
  return(pl)
}

SubEntities = c("Treg_out","IL10_out", "BBB",
                "Teff_in","IL17_in","IFNg_in","ODC irr. damaged")

Entities = c("Teff_out","Treg_out","Antigen" ,
             "EffectorMemory","NK_out","BBB",
             "IL10_out","IL17_out","IFNg_out", 
             "Treg_in", "Teff_in","ODC irr. damaged" ,
             "IL17_in","IFNg_in","IL10_in")

pl1<-Plotgeneration(DataFrame = DatiPlot, Stoch = TRUE, Mean=Mean, Entities=SubEntities)

ggsave(pl1,filename = paste0("Plot/StochDynamics.pdf"),device = "pdf",width = 10, height = 20 )


pl<-Plotgeneration(DataFrame = DatiPlot_Det, Stoch = FALSE, Mean=Mean, Entities=Entities)
ggsave(pl,filename = paste0("Plot/DetermDynamicsALL.pdf"),device = "pdf",width = 15, height = 15 )

pl<-Plotgeneration(DataFrame = DatiPlot, Stoch = TRUE, Mean=Mean, Entities=Entities)
ggsave(pl,filename = paste0("Plot/StochDynamicsALL.pdf"),device = "pdf",width = 10, height = 20 )
