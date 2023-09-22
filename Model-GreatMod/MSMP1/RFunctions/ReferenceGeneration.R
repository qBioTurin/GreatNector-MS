library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
## reading reference  generated from Datainvestigation.Rmd:
reference <- read_table2("~/temporalanalysis_MS/Model/Input/reference.csv")

## creating data.frame from literature:
ATZreferenceC1 = data.frame(Time = c(0,5,10,15,20,25,30,30*3,30*6) + 6*30,
														ATZ = c(0, 2800, 1250,450,200,150,100,50,0) )
ATZreferenceC2 = data.frame(Time = c(0,3,5,7,10,30,30*6) + 18*30,
														ATZ = c(0,2100,1600,400,100,50,0) )
ATZreference = rbind(ATZreferenceC1,ATZreferenceC2)
ATZreference$ATZ = ATZreference$ATZ /10^9/ 145453.8 * 6.022e+23

# CD3reference =  data.frame(Time = c(0,1,3,6,9,12,13,15,18,21,24,48,60,72,84) ,
# 													 CD3 = c(2,0.2,0.7,0.8,1,1.2,0.35,0.5,0.8,1,1.2,rep(2,4))*10^3 )

for(k in unique(reference$ClusterKmeans)){
	referenceG = reference %>%
		filter(ClusterKmeans == k) %>%
		dplyr::select(-ClusterKmeans) %>%
		group_by(ID) %>%
		dplyr::mutate(Time = Time * 30) %>%
		spread(Exp, Observation) %>%
	  dplyr::select(-CD4) %>%
		gather(Exp,Obs, Th1,Th17,Treg) %>%
		arrange(Time) %>%
	  ungroup() %>%
	  na.omit()
	
print(
	ggplot() +
		geom_line(data = referenceG, aes(x = Time, y = Obs, group = ID,col=ID)) +
		geom_point(data = referenceG, aes(x = Time, y = Obs, group = ID,col=ID)) +
		geom_line(data = referenceG %>%
								group_by(Exp,Time) %>%
								dplyr::summarise(Mean = mean(Obs)),
							aes(x = Time, y = Mean, linetype = Exp)) +
		facet_wrap(~Exp,scales = "free")
)
	Ref = rbind(referenceG ,
				 data.frame(Time = ATZreference$Time, ID = "ATZ", Exp= "ATZ", Obs = ATZreference$ATZ))
	RefMean = rbind(referenceG 	%>%
										group_by(Exp,Time) %>%
										dplyr::summarise(Mean = mean(Obs)) %>% as.data.frame(),
							data.frame(Exp= "ATZ", Time = ATZreference$Time, Mean = ATZreference$ATZ,stringsAsFactors = F))

	write.table(Ref,
							file = paste0("~/temporalanalysis_MS/Model/Input/referenceGroup",k,".csv"),
							quote = F, row.names = F,col.names = F)

	write.table(t(RefMean),
							file = paste0("~/temporalanalysis_MS/Model/Input/referenceMeanGroup",k,".csv"),
							quote = F, row.names = F,col.names = F)

}
