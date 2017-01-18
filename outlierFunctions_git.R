
# Loading necessary R packages
require(rrcov); require(MASS); require(reshape);require(ggplot2)

#####################################################################################################
## function to calculate cov^(-1/2) ##
#####################################################################################################
"%^%" <- function(x, n){
  with(eigen(x), vectors %*% (values^n * t(vectors)))
}

#####################################################################################################
## Calculating the Robust MCD-based Mahalanobis Distance Outlyingness measure ##
#####################################################################################################
getRMDO <- function(
  dat,
  nsamp="best", # the nsamp variable for CovMcd
  myalpha = 0.5 # alpha for proportion of obs to include in MCD-MD
){
  # Using fast-MCD to calculate robust mean and covariance
  mcd1 <- CovMcd(dat,alpha=myalpha,nsamp=nsamp) 
  x1 <- mcd1@center
  cn <- mcd1@cov
  
  # Calculating the Mahalanobis distance (MD), and then robust MD outlyingness (RMDO) based on robust mean and covariance
  md <- apply(dat,1,function(i){t(i-x1)%*%solve(cn)%*%(i-x1)})
  rmde <- 1 + sqrt(md)
  rmde1 <- 1/rmde
  rmdo <- 1-rmde1 # outlyingness values
#   
#   used <- rownames(dat)[mcd1@best] # the subjects who were included in the MCD calculations
#   
  # Based on 2016 paper: Contributions to Quadratic Form
  smpsd <- sqrt(diag(cn))
  d <- diag(1/smpsd)
  dsd <- d %*% cn %*% d
  
  trans1 <- sweep(dat,2,x1)
  corrmax <- dsd %^% (-1/2) %*% d
  contrib <- apply(trans1,1,function(i){corrmax %*% as.matrix(i)})
  w2 <- t(contrib**2)
  colnames(w2) <- colnames(dat)
  
  return(list(RMDO=rmdo, mcdCenter=x1, mcdCov=cn, initDat=dat, myalpha=myalpha, W2=w2))
  # mcdUsed=used, DSD=dsd, CorrMax=corrmax, 
}



#####################################################################################################
## Identifying Cut-off for Outliers with simulation ##
#####################################################################################################
idOutliers <- function(
  obsOV, # output from the getRMDO function
  rmdo_nsamp = 500, # number of iterations - same as nsamp in getRMDO.  "best" would be long but potentially do more than default of 500
  nsamps=250, # number of elements within each sample
  nreps=100, # number of samples to take and average over
  thisquant=0.95 # the quantile to use as cut-offs for outliers
){
  # Pulling necessary output from getRMDO
  mcdMn <- obsOV[["mcdCenter"]]
  mcdCov <- obsOV[["mcdCov"]]
  myalpha <- obsOV[["myalpha"]]
  
  # simulating multivariate normal distributions with the robust mean and covaraince from getRMDO
  outquants <- c()
  for(i in 1:nreps){
    tsamp<- mvrnorm(nsamps,mcdMn,mcdCov)
    trmdo <- getRMDO(tsamp,nsamp=rmdo_nsamp,myalpha=myalpha)   # this is where time is problematic
    outquants[i] <- quantile(trmdo[["RMDO"]],thisquant) # saving the quantile value only
  }
  
  # pulling the mean of all simulated quantiles -- will be the cut-off
  simquant <- mean(outquants)

  return(list(cutoff=simquant,cutoff_all=outquants))
}

#####################################################################################################
## Identifying Multivariate Outliers Function ##
#####################################################################################################

multiOut <- function(
  dat, # the data frame, excluding identifiers
  rmdo_alpha = 0.8, # see getRMDO -- myalpha
  rmdo_nsamp = "best",
  nsims = 100, # then number of total simuations (so nsims/iters per iteration)
  nsamps = 250, # the number of points to simulate in process of determining cut-off for outlier
  alpha=0.01 # error to allow in CIs for quantile cut-offs
){
  
  # Giving column and row names if none are previously present
  if(is.null(colnames(dat))){
    colnames(dat) <- paste("Col",1:ncol(dat),sep="")
  } 
  if(is.null(rownames(dat))){
    rownames(dat) <- paste("Row",1:nrow(dat),sep="")
  } 
  
  # Calling getRMDO to calculate outlyingness
  rmdos <- getRMDO(dat,nsamp=rmdo_nsamp,myalpha=rmdo_alpha)
  
  # Output of getRMDO
  rmdoCalc <- rmdos[["RMDO"]]
  centers <- rmdos[["mcdCenter"]]
  scatters <- data.frame(rmdos[["mcdCov"]])
  contVar <- data.frame(rmdos[["W2"]])
  contVar$ID <- rownames(contVar)
  
  # Simulating to identify "normal range"
  idOut <- idOutliers(obsOV=rmdos,nsamps=nsamps,nreps=nsims,thisquant = 1-alpha)
  outCut <- idOut[["cutoff"]]

  # Reformatting contVar to long form for ease of interpretation
  contVar1 <- suppressMessages(melt(contVar))
  
  # Making final lists of outliers
  outs <- data.frame(ID=rownames(dat),RMDO=rmdoCalc,Outlier=ifelse(rmdoCalc >= outCut,1,0))
  outs1 <- outs[order(outs$RMDO,decreasing = T),]
  
  return(list(initDat=as.data.frame(dat),RMDO=outs1, RMDOalpha=rmdo_alpha, Contribs=contVar1, CutOff=outCut))
} # end of multiOut function


#####################################################################################################
## Printing the results of multiOut ##
#####################################################################################################

printOuts <- function(
  multiOutdat, # data output by multiOut function
  nCont = 3 # the number of top contributors to consider
){
  
  ######### Plot 1: Distribution of mean RMDOs
  outs1 <- multiOutdat[["RMDO"]]
  outCut <- multiOutdat[["CutOff"]]
  # Plot of decreasing RMDO values, with the cut-off identified
  outs1$ID <- factor(outs1$ID,level=outs1$ID)
  outs1$Outlier <- factor(outs1$Outlier,level=unique(outs1$Outlier))
  
  print(ggplot(outs1,aes(x=RMDO,y=rep(0,nrow(outs1)),color=Outlier)) + geom_point(size=4) + geom_vline(xintercept=outCut,linetype="dashed",size=0.5) + 
          geom_text(data=data.frame(outCut),aes(label="Cut-Off",x=outCut,y=0,fontface="bold"),size=5.5,color="black",nudge_x=0.008,nudge_y=-0.005) + 
          labs(title="Robust MCD-based Mahalanobis Distance Outlyingness\nMultivariate Outlier Detection",x="RMDO",y="") + 
          scale_y_continuous(limits = c(-0.007, 0.007)) +
          theme(axis.ticks.y=element_blank(),axis.text.y=element_blank(),plot.title=element_text(face="bold"),
                legend.text=element_text(size=13),axis.text.x=element_text(colour="black"),axis.title.x=element_text(size=15),
                legend.title=element_text(size=15)))
  
  ######### Plot 2: Distribution of Contribution for each observation
  Contribs <- multiOutdat[["Contribs"]]
  
  contribs1 <- merge(Contribs,outs1)
  contribs2 <- contribs1[order(contribs1$RMDO),]
  contribs2$ID <- factor(contribs2$ID, levels=as.character(unique(contribs2$ID)))

  outline <- nrow(outs1)-length(which(outs1$Outlier == 1))+0.5; outlab <- max(contribs2$value)
  print(ggplot(contribs2,aes(x=value,y=ID,col=variable)) + geom_point() + 
          labs(title="Variation in Contribution to RMDO",x="Contribution Value",y="Participant") + 
          geom_hline(yintercept = outline,linetype="dashed") + 
          geom_text(data=data.frame(outline),aes(label="Outliers",y=outline,x=outlab,fontface="bold"),color="black",nudge_x=-5,nudge_y=1) +
          theme(legend.position="none",plot.title=element_text(face="bold"),
                axis.text.y=element_text(colour="black"),axis.text.y=element_text(size=7.5)))
  
  ######### Plot 3: Distribution of Contribution for each observation - Anomalies only
  if(any(outs1$Outlier == 1)){
    contribs3 <- subset(contribs2,Outlier == 1)
    contribs3$ID <- factor(contribs3$ID, levels=as.character(unique(contribs3$ID)))
    
    print(ggplot(contribs3,aes(x=value,y=ID,col=variable)) + geom_point(size=3) + 
            labs(title="Variation in Contribution to RMDO - Outliers Only",x="Contribution Value",y="Participant") + 
            scale_x_continuous(limits=c(0,max(contribs3$value))) + 
            theme(legend.position="none",plot.title=element_text(face="bold"),axis.text.y=element_text(colour="black")))
    
  } # end of "yes outliers exist"
  
  ######### Plot 4+: 2x2 plots showing outliers
  cont1 <- split(Contribs,Contribs$ID)
  cont2 <- lapply(cont1,function(i){i[order(i[,"value"],decreasing = T),][1:3,]})
  topCont <- do.call(rbind,cont2)

  topCont1 <- merge(outs1,topCont)
  
  if(any(topCont1$Outlier == 1)){ # some data sets will not have outliers
    topCont2 <- topCont1[which(topCont1$Outlier == 1),]
    topCont2 <- topCont2[order(topCont2$RMDO,topCont2$value,decreasing = T),]
    colnames(topCont2)[which(colnames(topCont1)=="value")] <- "contrib"
    
    initdat <- multiOutdat[["initDat"]]
    
    # Adding observed values to the table
    
    topCont2$value <- as.matrix(apply(topCont2[,c("ID","variable")],1,function(i){initdat[which(rownames(initdat) == i[1]),which(colnames(initdat) == i[2])]}))
    tsub <- as.character(unique(topCont2$ID))
    
    # loop through anomalies to plot based on top contributors
    for(i in tsub){
      toptmp <- topCont2[which(topCont2$ID == i),]; toptmp$rank <- 1:nrow(toptmp)
      subtmp <- i; rmdotmp <- toptmp$RMDO[1]
      vartmp <- as.character(unique(toptmp$variable))
      
      sizetmp <- ifelse(rownames(initdat) == subtmp,5.5,4)
      shapetmp <- ifelse(rownames(initdat) == subtmp,17,16)
      
      if(nCont == 2){
        coltmp <- ifelse(rownames(initdat) == subtmp,"red","blue")
        print(ggplot(initdat,aes(x=initdat[,vartmp[1]],y=initdat[,vartmp[2]])) + geom_point(col=coltmp,size=sizetmp,shape=shapetmp) +
                labs(title=i,x=vartmp[1],y=vartmp[2]) + 
                theme(plot.title=element_text(face="bold"),axis.text=element_text(colour="black"),
                      axis.title=element_text(size=15)))
      }else if(nCont == 3){
        print(ggplot(initdat,aes(x=initdat[,vartmp[1]],y=initdat[,vartmp[2]],col=initdat[,vartmp[3]])) + geom_point(size=sizetmp,shape=shapetmp) +
                labs(title=i,x=vartmp[1],y=vartmp[2],col=vartmp[3]) + 
                theme(plot.title=element_text(face="bold"),axis.text=element_text(colour="black"),
                      axis.title=element_text(size=15),legend.title=element_text(size=15)))
      }else if(nCont == 4){
        print(ggplot(initdat,aes(x=initdat[,vartmp[1]],y=initdat[,vartmp[2]],col=initdat[,vartmp[3]])) + geom_point(size=sizetmp,shape=shapetmp) +
                labs(title=paste(i, "(1,2,3)",sep=" "),x=vartmp[1],y=vartmp[2],col=vartmp[3]) + 
                theme(plot.title=element_text(face="bold"),axis.text=element_text(colour="black"),
                      axis.title=element_text(size=15),legend.title=element_text(size=15)))
        
        print(ggplot(initdat,aes(x=initdat[,vartmp[1]],y=initdat[,vartmp[2]],col=initdat[,vartmp[4]])) + geom_point(size=sizetmp,shape=shapetmp) +
                labs(title=paste(i, "(1,2,4)",sep=" "),x=vartmp[1],y=vartmp[2],col=vartmp[4]) + 
                theme(plot.title=element_text(face="bold"),axis.text=element_text(colour="black"),
                      axis.title=element_text(size=15),legend.title=element_text(size=15)))
      }else{
        cat("Chosen nCont not set up\n")
        browser()
      }
    } # end of loop through anomalies
  } # end of "yes outliers exist
  return(topCont2)
} # end of printOuts

