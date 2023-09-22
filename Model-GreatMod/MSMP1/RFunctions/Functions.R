init_m <- function(n_file, optim_v=NULL, referencePath)
{
  
  reference<-read.table(referencePath)
  # reference<-read.table("./Input/referenceGroup4.csv")
  colnames(reference) = c("Time","ID","Place","Obs")
  reference[which(reference$Place == "Th1"),"Place"]  = "Teff_out_t1"
  reference[which(reference$Place == "Th17"),"Place"]  = "Teff_out_t17"
  reference[which(reference$Place == "Treg"),"Place"]  = "Treg_out"
  
  yini.names <- readRDS(n_file)
  
  yini <- rep(0,length(yini.names))
  dim(yini)<- c(1,length(yini.names))
  yini <- as.data.frame(yini)
  names(yini)<-yini.names
  
  # From Simona/Alessandro 's file:
  # group 4:
  yini["Teff_out_t1"] = mean(reference[which(reference$Place == "Teff_out_t1" & reference$Time ==  0),"Obs"])
  yini["Teff_out_t17"] = mean(reference[which(reference$Place == "Teff_out_t17" & reference$Time ==  0),"Obs"])
  yini["Treg_out"] = mean(reference[which(reference$Place == "Treg_out" & reference$Time ==  0),"Obs"])
  
  return(matrix(as.integer(yini), ncol = 1))
}

MappingParams<- function(name,params_name, optim_v=NULL )
{
  p.names <- readRDS(params_name)
  names(optim_v) = p.names
  
  p=optim_v[name]
  
  #p = runif(1,min = p*0.8,max = p*1.5)
  
  return(matrix(p, ncol = 1))
}


EventFun = function(marking,time){
  newmarking <- marking
  names(newmarking) <- readRDS("/home/docker/data/Input/NAMES.RDS")
  
  if(time == 180 ){
    newmarking["ATZ"] = marking["ATZ"] + 1.159241e+13
  }else if(time == 540 ){
    newmarking["ATZ"] = marking["ATZ"] + 8.694307e+12
  }
  
  return(newmarking)
}

msqd<-function(reference, output)
{
  
  if(length(table(diff(output$Time))) > 1){
    return(1e20)
  }
  if(max(output$Time) < 84*30){
    return(1e20)
  }
  #reference <- as.data.frame(t(read.csv("Input/referenceMeanGroup2.csv", header = FALSE, sep = "")))
  
  reference = data.frame(reference,stringsAsFactors = F)
  ### Reference:
  colnames(reference) = c("Place","Time","Obs")
  print(reference)
  reference$Place = as.character(reference$Place)
  reference$Time = as.numeric(as.character(reference$Time))
  reference$Obs = as.numeric(as.character(reference$Obs))
  reference[which(reference$Place == "Th1"),"Place"]  = "Teff_out_t1"
  reference[which(reference$Place == "Th17"),"Place"]  = "Teff_out_t17"
  reference[which(reference$Place == "Treg"),"Place"]  = "Treg_out"
  
  
  PlaceToConsider = unique(reference$Place)
  
  subOutput = output[output$Time %in% reference$Time,c("Time",PlaceToConsider)]
  subOutput2 = do.call("rbind",
                       lapply(colnames(subOutput)[-1], function(x)
                         data.frame(Place = x, Time = subOutput$Time, Mark = subOutput[,x] )) )
  
  RefOutput = merge(reference,
                    subOutput2,
                    by=c("Time","Place"))
  RefOutput$Error = abs(RefOutput$Obs - RefOutput$Mark)
  RefOutput$ErrorRel = RefOutput$Error
  RefOutput$ErrorRel[RefOutput$Obs>1 ] = RefOutput$Error[RefOutput$Obs>1 ]/RefOutput$Obs[RefOutput$Obs>1 ]
  RefOutput$ErrorRel[RefOutput$Obs<=1 ] = abs(RefOutput$Error[RefOutput$Obs<=1 ])
  
  ErrorRef = sum(do.call("rbind",
                         lapply(unique(RefOutput$Place), function(x)
                           mean(RefOutput[RefOutput$Place==x,"ErrorRel"])
                         )
  )
  )
  
  ## Antigen and damage at minimum
  MeanAnt = mean(output[,"Antigen"])
  
  # Error definition:
  Error = mean(2*ErrorRef + MeanAnt )
  
  return(Error)
}

#
# output<-read.table("./MS_Model_calibration/MS_Model-calibration-1.trace",header = TRUE)
# reference<-t(read.table("./Input/referenceMeanGroup3.csv"))
## OLD errors:
##Interl_out <- apply(output[,c("INFg_out","IL17_out","IL10_out")],2,"mean")
# Interl_out <- output[day.inf,c("INFg_out","IL17_out","IL10_out")]
# Interl_in <- output[day.inf+2,c("INFg_in","IL17_in","IL10_in")]
#
# # Out: number of INFg 42 IL17	8 IL10 13	0.59	0.76	1.03 (healthy)
# Error_Interl_out <- sum((Interl_out - reference[3:5] )^2 )
# # In: number of INFg 0.59 IL17	076 IL10 1.03 (healthy)
# Error_Interl_in <- sum((Interl_in - reference[6:8] )^2)
#
# # Number of Teff and Treg out: (I put the duoble weight!!)
# ErrorTcell <- (c(Teff_out[day.inf1],Treg_out[day.inf1] ) - reference[9:10])^2
# if(Teff_out[ntime]== 0) ErrorTcell = 10^6
# if(Treg_out[ntime]== 0) ErrorTcell = 10^6
#
# # and Teff should be twice the treg
# ErrorTcell2 = mean( (0.5* Teff_out[(day.inf + (0:5)*24)  ] - Treg_out[(day.inf + (0:5)*24)  ] )^2)
#
# # max value of Teff_out is at the 5th day
# if(!is.na(Teff_out[day.inf+5*24]))
# {
# 	Error = Error+ 2*(max(Teff_out) - Teff_out[day.inf+5*24])^2
# }
# if(!is.na(A[day.inf+7*24]))
# {
# 	Antigene_left_middle<-A[day.inf+(3)*24]
# 	Antigene_left<-mean(A[day.inf+(7:14)*24])
# 	Error = Error+ (max(A)*.95 - Antigene_left )^2 + (max(A)*.5 -Antigene_left_middle)^2
# 	# After 7 days I want that the antigen is 90% reduced
#
# }
# # I want to minimize the ODC damage
# # ODC_error= mean(output[,"ODC_le1"])^2
#
# # Error = Error +  perc_Teff + perc_Treg +  sum(Error_Interl_out) + sum(Error_Interl_in)
# #Error = Error +  perc_Teff + perc_Treg + sum(ErrorTcell) + ODC_error + ErrorTcell2 + Error_Interl_out
# Error = Error +  perc_Teff + perc_Treg + sum(ErrorTcell) + ErrorTcell2 + Error_Interl_out + Error_Interl_in

