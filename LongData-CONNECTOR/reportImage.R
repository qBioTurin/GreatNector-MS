library(connector)
library(ggplot2)
library(patchwork)
library(dplyr)
library(viridis)

# CONNECTORList = S.cl$CONNECTORList
# clusterData = S.cl
# g = cluster[i]

reportPlot.generation = function(CONNECTORList,
                                 clusterData,
                                 g ,
                                 title = NULL){
  
  ## Plot: data overview
  ### length curves distribution
  pl00.1 = ggplot(data.frame(l = CONNECTORList$LenCurv) ) +
    geom_bar(aes(x = l.Freq)) + 
    theme_bw() +
    geom_text(aes(x= min(CONNECTORList$LenCurv)+1,
                  y = max(table(CONNECTORList$LenCurv))+1),
              label = paste0("Total curves = ", length(CONNECTORList$LenCurv) )) +
    labs(x = "Number of points",
         y= "Number of curves " )
  
  ### Time point distribution
  
  pl00.2 = ggplot( CONNECTORList$Dataset ) +
    geom_bar(aes(x = Time)) + 
    theme_bw() +
    labs(x = "Time",
         y= "Number of curves" )
  
  ## Plot: data visualization
  pl0 = TimeGridDensity(CONNECTORList)$TimeGrid_plot
  
  ## Plot: CrossLogLikelihood
  CrossLogLike<-BasisDimension.Choice(CONNECTORList,2:8,Cores = 10)
  pl1 = CrossLogLike$CrossLogLikePlot
  pl2 = CrossLogLike$KnotsPlot
 
    
    ## Plot: fDB indexes
    IndexesPlot.Extrapolation(clusterData,q = 0.0)-> indexes
  pl3.1 = ggplot(indexes$Plot$plot_env$Indexes.Rep %>% 
                   filter(Index == "fDB")) +
    facet_wrap(~Index, 
               scales = "free") + 
    geom_violin(aes(x = Cluster, y = V, 
                    fill = ClusterH, group = Cluster), scale = "width") + 
    geom_line(data = indexes$Plot$plot_env$Indexes.MostProb %>% filter(Index == "fDB"),
              aes(x = Cluster, y = V, col = "Most probable")) + 
    geom_jitter(aes(x = Cluster, 
                    y = V),
                color = "black", 
                width = 0.1, height = 0,
                alpha = 0.5) + 
    scale_fill_viridis("", discrete = TRUE, alpha = 0.6) + 
    labs(size = "Counts freq.", col = "", x = "Number of Clusters", 
         y = "Indexes Values") + 
    scale_x_continuous(breaks = unique(indexes$Plot$plot_env$Indexes.MostProb$Cluster)) + 
    scale_color_manual(values = c(`Most probable` = "blue")) + 
    theme(axis.text = element_text(size = 14, hjust = 0.5), 
          axis.text.x = element_text(vjust = 0.5, hjust = 1), 
          axis.title = element_text(size = 16, face = "bold"), 
          axis.line = element_line(colour = "black"),
          plot.title = element_text(size = 30, 
                                    face = "bold", 
                                    vjust = 1, 
                                    lineheight = 0.6), 
          legend.text = element_text(size = 16),
          legend.position = "bottom", 
          legend.key = element_blank(),
          legend.title = element_text(size = 16, 
                                      face = "bold"), 
          legend.key.size = unit(0.9, "cm"), 
          legend.key.width = unit(0.9, "cm"),
          panel.background = element_rect(colour = NA), 
          plot.background = element_rect(colour = NA),
          plot.margin = unit(c(5, 
                               5, 5, 5), "mm"),
          strip.text = element_text(size = 20))+ 
    coord_cartesian(ylim=c(0, 5))
  
  pl3.2 = ggplot(indexes$Plot$plot_env$Indexes.Rep %>% filter(Index == "Tightness")) +
    facet_wrap(~Index, 
               scales = "free") + 
    geom_violin(aes(x = Cluster, y = V, 
                    fill = ClusterH, group = Cluster), scale = "width") + 
    geom_line(data = indexes$Plot$plot_env$Indexes.MostProb %>% filter(Index == "Tightness"),
              aes(x = Cluster, y = V, col = "Most probable")) + 
    geom_jitter(aes(x = Cluster, 
                    y = V),
                color = "black", 
                width = 0.1, height = 0,
                alpha = 0.5) + 
    scale_fill_viridis("", discrete = TRUE, alpha = 0.6) + 
    labs(size = "Counts freq.", col = "", x = "Number of Clusters", 
         y = "Indexes Values") + 
    scale_x_continuous(breaks = unique(indexes$Plot$plot_env$Indexes.MostProb$Cluster)) + 
    scale_color_manual(values = c(`Most probable` = "blue")) + 
    theme(axis.text = element_text(size = 14, hjust = 0.5), 
          axis.text.x = element_text(vjust = 0.5, hjust = 1), 
          axis.title = element_text(size = 16, face = "bold"), 
          axis.line = element_line(colour = "black"),
          plot.title = element_text(size = 30, 
                                    face = "bold", 
                                    vjust = 1, 
                                    lineheight = 0.6), 
          legend.text = element_text(size = 16),
          legend.position = "bottom", 
          legend.key = element_blank(),
          legend.title = element_text(size = 16, 
                                      face = "bold"), 
          legend.key.size = unit(0.9, "cm"), 
          legend.key.width = unit(0.9, "cm"),
          panel.background = element_rect(colour = NA), 
          plot.background = element_rect(colour = NA),
          plot.margin = unit(c(5, 
                               5, 5, 5), "mm"),
          strip.text = element_text(size = 20))
  
  pl3 =pl3.1 + pl3.2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  ## Plot:   Consensus matrix
  ConsMatrix.Extrapolation(clusterData)-> ConsInfo
  pl4 = ConsInfo[[paste0("G",g)]]$ConsensusPlot
  
  ## Plot: 
  MostProbableClustering.Extrapolation(clusterData,g,q = .0) -> MostProbableClustering
  
  Cl = MostProbableClustering$FCM$cluster$ClustCurve
  Cl$Cluster = MostProbableClustering$FCM$cluster$cluster.names[Cl$Cluster]
  meancurves = MostProbableClustering$FCM$prediction$meancurves %>% as.data.frame()
  colnames(meancurves) = MostProbableClustering$FCM$cluster$cluster.names
  
  meancurves$Time = MostProbableClustering$FCM$TimeGrid
  meancurves = meancurves %>% tidyr::gather(-Time,key = "Cluster",value = "Mean")
  
  pl5 = ggplot() +
    geom_line(data = Cl,aes(x = Time, y = Observation, group = ID, col = Cluster) ) +
    geom_line(data = meancurves,aes(x = Time, y = Mean) ) +
    facet_wrap(~Cluster,nrow = 2)+
    theme_bw()
  
  pl00 = (pl00.1|pl00.2|pl0) +
    plot_annotation(title = 'A) Data overview',
                    theme = theme(plot.title = element_text(size = 18,face = "bold"))) 
  
  pl01 = (pl2 + pl1) + 
    plot_layout(widths = c(2, 1))+
    plot_annotation(title = 'B) Selection of p',
                    theme = theme(plot.title = element_text(size = 18,face = "bold"))) 
  
  pl02 = pl3 +
    plot_annotation(title = 'C) Selection of G',
                    theme = theme(plot.title = element_text(size = 18,face = "bold")))
  
  pl03 = (pl4 | pl5)+
    plot_annotation(title = 'D) Clustering',
                    theme = theme(plot.title = element_text(size = 18,face = "bold")))
  
  pl = (  wrap_elements(pl00) /  wrap_elements(pl01) /  wrap_elements(pl02) /  wrap_elements(pl03) )+ 
    plot_layout(nrow = 4, heights = c(1, 2, 1, 1))
  
  return(pl)
}

################################
load("RData/Feature.RData")

exper = c("CD4","Th17_prod","Th1_prod","Treg_prod")
cluster = c(3,4,3,3,3)
names(cluster) = exper

for( i in exper){
  if(i == "Th1_prod")
    load( paste0("RData/Cluster",i,"_p4.RData") )
  else
    load( paste0("RData/Cluster",i,"_p5.RData") )
  
  
  pl = reportPlot.generation(S.cl$CONNECTORList,S.cl,cluster[i]) Tre+
        plot_annotation(title = paste0("CONNECTOR analysis: ",gsub(i,pattern = "_prod",replacement = "")),
                  theme = theme(plot.title = element_text(size = 20,face = "bold")))
  
  ggsave(pl,filename = paste0("Plot/SupplementaryFigure",gsub(i,pattern = "_prod",replacement = ""),"Indexes.pdf"),
         device = "pdf",height = 25,width = 15)
  
}



