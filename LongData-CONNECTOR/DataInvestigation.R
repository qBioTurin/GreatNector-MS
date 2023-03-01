library(connector)
library(dplyr)
library(tidyverse)
library(ggplot2)

#exper = c("CD4","EDSS","Th17","Th1","Treg"),"Th17Treg"
exper = c("CD4","EDSS","Th17_prod","Th1_prod","Treg_prod")

load("RData/Feature.RData")
cluster = c(3,4,3,3,3)
names(cluster) = exper
ClList = list()
TempList = list()
meancurvesList= list()
for( i in exper){
	if(i == "Th1_prod")
		load( paste0("RData/Cluster",i,"_p4.RData") )
	else
		load( paste0("RData/Cluster",i,"_p5.RData") )

	load( paste0("RData/temporal",i,".RData") )
	print(i)
	### cluster
	MostProbableClustering.Extrapolation(S.cl,cluster[i]) -> MostProbableClustering
	Cl = MostProbableClustering$FCM$cluster$ClustCurve %>% select(ID,Cluster) %>% distinct()
	Cl$Cluster = MostProbableClustering$FCM$cluster$cluster.names[Cl$Cluster]
	meancurves = MostProbableClustering$FCM$prediction$meancurves %>% as.data.frame()
	colnames(meancurves) = MostProbableClustering$FCM$cluster$cluster.names

	meancurves$Time = MostProbableClustering$FCM$TimeGrid
	meancurves = meancurves %>% tidyr::gather(-Time,key = "Cluster",value = "Mean")
	meancurves$Exp = i

	colnames(Cl) = c("ID",paste0("Cluster",i))
	ClList[[i]] = Cl

	### temporal
	colnames(Cl) = c("ID","Cluster")
	Cl = merge(MostProbableClustering$CONNECTORList$LabCurv,Cl)
	Cl = Cl %>% rename( IDconnector = ID , ID = IDold)

	temporalCurve = merge(temporalCurve,
												Cl) # %>% select(ID,IDconnector,Cluster) )
	temporalCurve$Exp = i
	TempList[[i]] = temporalCurve
	meancurvesList[[i]] = meancurves
}

meancurves = do.call("rbind",meancurvesList)
saveRDS(meancurves, file = "RData/MeanCurves.Rds")

TemporalCurve = do.call("rbind",TempList)
TemporalCurve$Observation = as.numeric(TemporalCurve$Observation )
ID = TemporalCurve %>% select(ID, IDconnector) %>% distinct()

Clusters = ClList %>% reduce(full_join, by = "ID" )
Clusters = Clusters %>% rename( IDconnector = ID )
Clusters = merge(ID,Clusters)
CompleteFeat = merge(Feature,Clusters)
TemporalCurve$ID = as.factor(TemporalCurve$ID)

rel = grep(pattern = "RELAPSES", colnames(TemporalCurve))
TemporalCurve$NumberRelapses = 0

for( i in unique(TemporalCurve$ID)){
	TemporalCurve[which(TemporalCurve$ID == i),"NumberRelapses"] =
		length(TemporalCurve[which(TemporalCurve$ID == i)[1],rel][!is.na(TemporalCurve[which(TemporalCurve$ID == i)[1],rel] )]
												 )
}

TemporalCurve = TemporalCurve %>%
	mutate(Relapses = ifelse(NumberRelapses ==0, "No rel","At least one"))
TemporalCurve$TypeRel = sapply(1:length(TemporalCurve[,1]),
															 function(x) gsub(pattern = "NA",replacement = "",paste(TemporalCurve[x,rel],collapse = "")) )

TemporalCurve$NumberRelapses = as.factor(TemporalCurve$NumberRelapses)

saveRDS(TemporalCurve, file = "RData/Summary.Rds")

ggplot(TemporalCurve%>%
			 	arrange(Time) %>% drop_na(Observation) ) +
	geom_point(aes(x = Time, y = Observation, group = ID, col = Relapses, shape = ID))+
	scale_shape_manual(values=1:nlevels(TemporalCurve$ID)) +
	geom_line(aes(x = Time, y = Observation,group = ID, col = Relapses))+
	facet_grid(TypeRel ~ Cluster) +
	theme(legend.position = "bottom") +
	theme_bw()


