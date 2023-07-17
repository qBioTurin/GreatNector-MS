library(dplyr)
library(readxl)
library(tidyr)
library(connector)

ALEMTUZUMAB_DATA <- read_excel("Data/ALEMTUZUMAB_DATA_Connector.xlsx", col_types = c())
length( unique( ALEMTUZUMAB_DATA$ID ) )

unique(ALEMTUZUMAB_DATA_OLD$`experimental code`)[ unique( ALEMTUZUMAB_DATA_OLD$`experimental code` ) %in% unique( ALEMTUZUMAB_DATA$ID ) ]


#### Feature file generation

Feature = ALEMTUZUMAB_DATA %>% dplyr::select_if(!grepl("EDSS", names(.)))
save(Feature,file = "RData/Feature.RData")

#### Time Series generation
## EDSS

temporalCurve = ALEMTUZUMAB_DATA %>%
  dplyr::select_if(grepl("ID|EDSS", names(.))) %>%
	gather("Step","Observation",-ID)

temporalCurve$Step = gsub(temporalCurve$Step,
                          pattern = "(EDSS AT )|(EDSS)",
                          replacement = "")

steps = unique(temporalCurve$Step)
time = 0:( length(steps) - 1)
names(time) = steps

temporalCurve$Time = time[temporalCurve$Step]
temporalCurve = temporalCurve %>% 
  dplyr::select(ID,Observation,Time) %>%
  dplyr::mutate( Time = ifelse(Time < 5, Time*6, (Time-2)*12 ) )

save(temporalCurve,file = "RData/temporalEDSS.RData")
CD3reference =  data.frame(Time = c(0,1,3,6,9,12,13,15,18,21,24,36,48,60,72) ,
													 CD3 = c(2,0.2,0.7,0.8,1,1.2,0.35,0.5,0.8,1,1.2,rep(2,4))*10^3 )

## Others

for( i in c("CD4","Treg","Th1","Th17")){
	ALEMTUZUMAB_DATA <- read_excel("Data/ALEMTUZUMAB DATA Connector_Puliti.xlsx",
																 sheet = i)

	colnames(ALEMTUZUMAB_DATA) = paste0("M",colnames(ALEMTUZUMAB_DATA))
	temporalCurve = ALEMTUZUMAB_DATA %>%
	  rename( ID = Mmonths  ) %>%
		gather("Step","Observation",-ID)

	steps = unique(temporalCurve$Step)
	time = 0:(length(steps) - 1)
	names(time) = steps

	temporalCurve$Time = time[temporalCurve$Step]
	temporalCurve = temporalCurve %>% dplyr::select(ID,Observation,Time)%>%
		dplyr::mutate(Time = ifelse(Time < 5, Time*6, (Time-2)*12 ),
		              Observation = as.numeric(Observation) )

	if(i != "CD4"){
		temporalCurve -> temporalCurve.subPerc
		load("RData/temporalCD4.RData")
		colnames(temporalCurve) = c("ID","ObservationCD4","Time")
		temporalCurve <- merge(temporalCurve,temporalCurve.subPerc,by = c("ID","Time"))
		temporalCurve = temporalCurve %>% drop_na()
		temporalCurve$ObservationCD4 = as.numeric(temporalCurve$ObservationCD4)
		temporalCurve$Observation = as.numeric(temporalCurve$Observation)
		temporalCurveTMP = temporalCurve
		temporalCurve$prop = temporalCurve$Observation/temporalCurve$ObservationCD4 *100
		temporalCurve = temporalCurve %>% dplyr::select(ID,prop,Time) %>% rename(Observation = prop)
		#save(temporalCurve,file = paste0("RData/temporal",i,"_prop.RData") )

		temporalCurveTMP -> temporalCurve
		temporalCurve$prod = temporalCurve$Observation * temporalCurve$ObservationCD4 /100
		temporalCurve = temporalCurve %>% dplyr::select(ID,prod,Time) %>% rename(Observation = prod)
		save(temporalCurve,file = paste0("RData/temporal",i,"_prod.RData") )

	}else{
		temporalCurve2 = merge(CD3reference,temporalCurve, by = "Time",all.y = T)
		temporalCurve = temporalCurve2 %>%
			dplyr::mutate(Observation = as.numeric(Observation)/100 * CD3) %>%
			dplyr::select(-CD3) %>%
			dplyr::select(ID,Observation,Time)
		save(temporalCurve,file = paste0("RData/temporal",i,".RData") )
	}

}

#############################################
### Th17/Treg
#############################################

load("./RData/temporalTh17_prod.RData")
colnames(temporalCurve) = c("ID","ObservationTh17","Time")
temporalCurve -> temporalCurveTh17

load("RData/temporalTreg_prod.RData")
colnames(temporalCurve) = c("ID","ObservationTreg","Time")
temporalCurve -> temporalCurveTreg

temporalCurve = merge(temporalCurveTreg,temporalCurveTh17,by = c("ID","Time"))
temporalCurve = temporalCurve %>% drop_na() 
temporalCurve$ObservationTreg = as.numeric(temporalCurve$ObservationTreg)
temporalCurve$ObservationTh17 = as.numeric(temporalCurve$ObservationTh17)
temporalCurve$Observation = temporalCurve$ObservationTh17/temporalCurve$ObservationTreg
temporalCurve = temporalCurve %>% dplyr::select(ID,Observation,Time)
save(temporalCurve,file = paste0("RData/temporalTh17Treg.RData") )
