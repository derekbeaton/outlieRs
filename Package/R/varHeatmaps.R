varHeatmaps <- function(
  dat # output from multiOut$hidden_detail (or any list with contrib, and CorrMat elements)
){
  contrib <- dat$Contribs
  corrmat <- dat$corrMat
  
  # heatmap of outliers by variable contribution -- hoping to see which variables contribute the most to each outlier
  contribOut <- contrib[which(contrib$ID %in% as.character(dat$RMDO$ID[which(dat$RMDO$Outlier == 1)])),]
  print(ggplot(contribOut,aes(x=variable,y=ID,fill=value)) + geom_raster() + scale_fill_gradientn(colours = heat.colors(100)) + 
          labs(title=sprintf("Heatmap of Variable Contribution to Outlier"),fill="Score") + 
          theme(plot.title=element_text(face="bold"),axis.text.x = element_text(angle=90,hjust=1),
                axis.title=element_blank()))
  
  # heatmap variable by variable corrmax transformation  -- can see if covariance of some variables dominate
  corrmat1 <- suppressMessages(melt(data.frame(v1=rownames(corrmat),corrmat)))
  colnames(corrmat1)[2] <- "v2"
  corrmat1$v1 <- factor(as.character(corrmat1$v1), levels = unique(corrmat1$v1[order(as.character(corrmat1$v1))]))
  corrmat1$v2 <- factor(as.character(corrmat1$v2), levels = unique(corrmat1$v2[order(as.character(corrmat1$v2))]))
  
  print(ggplot(corrmat1,aes(x=v1,y=v2,fill=value)) + geom_raster() + scale_fill_gradientn(colours = heat.colors(100)) + 
          labs(title=sprintf("Heatmap of CorrMax Transformation"),fill="Score") + 
          theme(plot.title=element_text(face="bold"),axis.text.x = element_text(angle=90,hjust=1),
                axis.title=element_blank()))
}
