library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

plot.generation = function(folder, sensFolder = NULL, reference = NULL, All = F,event_Antig,event_ATZ){

  filesAll <- list.files(paste0("./",folder,"/"),pattern='*.trace')
  filesTOread = filesAll[1:(min(length(filesAll),100))]
  
  ## Calibration input files:
  CalibF = F
  Calib_fileOptim <- list.files(paste0("./",folder,"/"),pattern='*calibration.RData',full.names = T)
  Calib_fileConfig <- list.files(paste0("./",folder,"/"),pattern='*-calibration_optim-config.csv',full.names = T)
  if(length(Calib_fileConfig)>0 & length(Calib_fileOptim)>0){
    load(Calib_fileOptim[1])
    calibration_optim_trace <-read.csv(Calib_fileConfig,sep = "")
    optim_trace = calibration_optim_trace[which.min(calibration_optim_trace[,"distance"]),"id"]
    optim_trace_file = grep(pattern = paste0("*-",optim_trace,".trace"),x = filesAll,value = T)
    OptimCalib =	data.frame(ID = optim_trace, read.csv(paste0('./',folder,'/',optim_trace_file), sep="")) %>%
      gather(Place, Obs,-ID, -Time)
    CalibF = T
  }

  if(!is.null(sensFolder)){
    sens_fileOptim <- list.files(paste0("./",sensFolder,"/"),pattern='ranking_',full.names = T)
    load(sens_fileOptim)
    rank = rank[order(rank$measure),]
    bestID = as.numeric(gsub(".*?([0-9]+).*", "\\1",  rank$id[1])) 
    filesTOread = rank$id[1:(min(length(filesAll),100))]
  }
  
	traceAll <- lapply(filesTOread, function(x) {
		IDsim<- as.numeric(gsub(pattern="(.trace)", gsub(pattern="([[:graph:]]+(-){1})", x=x, replacement=""), replacement="") )
		data.frame(ID = IDsim, read.csv(paste0('./',folder,'/',x), sep="",stringsAsFactors = F))
	} )
	trace<-do.call("rbind", traceAll)
	traceDF =	trace %>% 
	  mutate(Teff_Treg = Teff_in/(Treg_in + 1e-6) )%>%
	  dplyr::select(-Teff_in,-Treg_in) %>%
	  gather(Place, Obs,-ID, -Time) 
	 

	if(length(table(diff(trace$Time))) > 1)
		warning("Lsoda does not converge!!!!!")

	if(All){

		Entities = c("Teff_out_t1","Teff_out_t17","Treg_out",
								 "Antigen" , "ATZ","Teff_Treg")
		traceDF =	traceDF %>%
			filter(Place %in% Entities)
		if(CalibF){
		  OptimCalib =	OptimCalib %>%
		    filter(Place %in% Entities)
		}
	}

	traceDF = traceDF %>% mutate(Time = Time - 30*6)
	
	if(CalibF){
	  pl1 = ggplot() 	+
	    #geom_line(data = traceDF, aes(x = Time, y = Obs, group = ID),col = "grey",alpha=0.2) +
	    geom_line(data = OptimCalib,
	              aes(x = Time- 30*6, y = Obs)) +
	    facet_wrap(~Place,scales = "free")
	}else if(!is.null(sensFolder)){
	  pl1 = ggplot() 	+
	    geom_line(data = traceDF, aes(x = Time, y = Obs, group = ID),col = "grey",alpha=0.2) +
	    geom_line(data = traceDF %>%
	                filter(ID == bestID),
	              aes(x = Time, y = Obs)) +
	    facet_wrap(~Place,scales = "free")
	}else{
	  pl1 = ggplot() 	+
	    geom_line(data = traceDF, aes(x = Time, y = Obs, group = ID),col = "grey") +
	    geom_line(data = traceDF %>%
	                group_by(Place,Time) %>%
	                dplyr::summarise(Mean = mean(Obs)),
	              aes(x = Time, y = Mean)) +
	    facet_wrap(~Place,scales = "free")
	}


	steps = seq(min(traceDF$Time),max(c(traceDF$Time)),by = 30*6)

	### reference

	if(!is.null(reference) ){
		referenceG <- as.data.frame(read.csv(reference,
																				header = FALSE,
																				sep = "",stringsAsFactors = F))
		colnames(referenceG) = c("Time","ID","Place","Obs")
		#referenceG[which(referenceG$Place != "ATZ"),"Time"] = referenceG[which(referenceG$Place != "ATZ"),"Time"]+ 180

		referenceG[which(referenceG$Place == "Th1"),"Place"]  = "Teff_out_t1"
		referenceG[which(referenceG$Place == "Th17"),"Place"]  = "Teff_out_t17"
		referenceG[which(referenceG$Place == "Treg"),"Place"]  = "Treg_out"

		referenceG = referenceG %>% filter(Time <= max(traceDF$Time)) %>% 
		  mutate(Time = ifelse(Place == "ATZ", Time - 6*30, Time ))
   
		pl1 = pl1 +
			geom_boxplot(data = referenceG, aes(x = Time, y = Obs, group = Time,col="Ref"),alpha = 0.2) +
			geom_point(data = referenceG, aes(x = Time, y = Obs, group = ID,col="Ref"),alpha = 0.2) 
	}

	events = data.frame(times = c(event_Antig,event_ATZ),
											type = c(rep("Antigen",length(event_Antig)),rep("ATZ",length(event_ATZ))) )
	pl1 = pl1 +
		scale_x_continuous(name = "Months",breaks = steps,labels = steps/30 ) +
		geom_vline(data= events, aes(col = type, xintercept = times),linetype = "dashed",alpha = .5)+
		scale_color_manual(values = c("Antigen" = "blue","ATZ" = "green","Ref" = "red")) +
		theme_bw()

return(list(pl = pl1,trace = trace))
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




