pfComposite=function(TR,
                     bins=NULL,
                     nboot=1000,
                     binning=TRUE,
                     conf=c(0.05,0.95))
{
  
  ## IF TR is a matrix
  if (is.matrix(TR) | is.data.frame(TR)){
    ID=unique(TR[,1])
    lgth=c()
    for (i in 1:length(ID)){
      lgth[i]=length(na.omit(TR[TR[,1]==ID[i],1]))
    }
    m=max(lgth)
    Age=matrix(nrow=m,ncol=length(ID))
    TransData=matrix(nrow=m,ncol=length(ID))
    for (i in 1:length(ID)) {
      Age[,i]=c(TR[TR[,1]==ID[i],3], rep(NA, m-length(TR[TR[,1]==ID[i],3])))
      TransData[,i]=c(TR[TR[,1]==ID[i],4], rep(NA, m-length(TR[TR[,1]==ID[i],4])))    
    } 
    colnames(TransData)=ID
    colnames(Age)=ID
    TR=structure(list(Age=structure(Age,class="matrix") ,TransData=structure
                      (TransData,class="matrix"),Method="unspecified"))
    
  }
  
  if(binning==TRUE){
    # Define the sequence for binning if unspecified
    if(is.null(bins)){
      AgeRes=matrix(nrow=length(TR$Age[,1])-1,ncol=length(TR$Age[1,]))
      for (i in 1:length(TR$Age[1,])){
        AgeRes[,i]=c(diff(TR$Age[,i]))
      }
      width=ceiling(median(na.omit(AgeRes))/10)*10
      binI=floor(min(na.omit(TR$Age))/10)*10
      binF=ceiling(max(na.omit(TR$Age))/10)*10
      bins=seq(binI,binF,width)
    }
    
    # If specified, values for plotting are:
    if(is.null(bins)==FALSE){
      width=bins[2]-bins[1]
      binI=bins[1]
      binF=bins[length(bins)]
    }
    
    # Matrix to store results
    result=matrix(ncol=length(TR$Age[1,]),nrow=length(bins)-1)
    #sm_result=matrix(ncol=length(IDn),nrow=length(bins)-1)
    
    ## Binning procedure
    for (k in 1:length(TR$Age[1,])){
      c1 <- cut(TR$Age[,k], breaks = bins)
      tmean=tapply(TR$TransData[,k], c1, mean)
      result[,k]=c(as.numeric(tmean))
    } 
    # suppres Inf values occuring with specific charcoal series (binary series)
    result[!is.finite(result)]=NA
  }
  
  ##
  if(binning==FALSE){
    bins=TR$Age[,1]
    binF=max(bins)
    binI=min(bins)
    width=bins[2]-bins[1]
    mboot=matrix(nrow=length(bins),ncol=nboot)
    result=TR$TransData
  } else mboot=matrix(nrow=length(bins)-1,ncol=nboot)
  
  ## Bootstrap procedure
  for (i in 1:nboot){
    ne=sample(seq(1,length(result[1,]),1),length(result[1,]),replace=TRUE)
    ## If has one bin
    if(dim(result)[1]==1){
      mboot[1,i]=mean(result[,ne],na.rm=TRUE)
    } else mboot[,i]=c(apply(result[,ne],1,mean,na.rm=TRUE))
  }
  bootci=t(apply(mboot, 1, quantile, probs = conf,  na.rm = TRUE))
  bootmean=t(apply(mboot, 1, mean,  na.rm = TRUE))
  
  # write out the transformed data
  
  if(binning==TRUE){centres=bins[1:length(bins)-1]+width/2} else centres=bins
  if(binning==FALSE){width=NA}
  
  result2=as.data.frame(cbind(centres,t(bootmean),bootci))
  colnames(result2)=c("AGE","MEAN",as.character(conf))
  colnames(result)=colnames(TR$TransData)
  
  output=structure(list(BinnedData=structure(result,row.names = as.character(centres),col.names=colnames(TR$TransData), class = "matrix"),
                        mboot=mboot,
                        BinCentres=centres,    
                        Result=result2,
                        BinWidth=width,
                        nboot=nboot,
                        binning=binning,
                        BootMean=structure(bootmean,row.names = as.character(centres),col.names=c("Mean"),class = "matrix"), 
                        BootCi=structure(bootci,row.names = as.character(centres),col.names=as.character(conf),class = "matrix") ))
  class(output)="pfComposite"
  return(output)  
}

######SUMMARY########

######PLOT########

plot.pfComposite=function(x,type="ci",conf=c(0.05,0.95),palette="jet",add="NONE",...){
  # Value for plotting:
  w=(x$BinCentres[2]-x$BinCentres[1])/2
  
  if (type=="ci"){
    
    if(add=="sitenum")
    par(mfrow=c(2,1))
   
    bootci1=t(apply(x$mboot, 1, quantile, probs = conf,  na.rm = TRUE))    
    
    plot(x$BinCentres,x$BootMean, xlim=c(max(x$BinCentres)+w,min(x$BinCentres)-w), ylim= c(min(bootci1,na.rm=T),max(bootci1,na.rm=T)), axes=F, mgp=c(2,0,0),
         main=paste("Composite"), font.main=1, lab=c(8,5,5), 
         ylab="Composite", xlab="Age (cal yr BP)", cex.lab=1, pch=16, cex=0.5, type="o")
    axis(1); axis(2, cex.axis=1)
    axis(side = 1, at = seq(0, 99000, by = 500), 
         labels = FALSE, tcl = -0.2)  
    for (i in 1:length(conf)){
      lines(x$BinCentres,bootci1[,i],lty=2)
      pos=which.min(is.na(bootci1[,i]))
      text(min(x$BinCentres)-200,bootci1[pos,i],paste(conf[i]*100,"%",sep=""),col="black")
    }
    
    # Plot site number
    if(add=="sitenum"){
    sitenum=length(x$BinnedData[1,])-rowSums(is.na(x$BinnedData))
    plot(x$BinCentres,sitenum,xlim=c(max(x$BinCentres)+w,min(x$BinCentres)-w), ylim= c(min(sitenum,na.rm=T),max(sitenum,na.rm=T)), axes=F, mgp=c(2,0,0),
         main=paste("Sites #"), font.main=1, lab=c(8,5,5), 
         ylab="Sites #", xlab="Age", cex.lab=0.8, pch=16, cex=0.5, type="o")
    axis(1); axis(2, cex.axis=1)
    axis(side = 1, at = seq(0, 99000, by = 500), 
         labels = FALSE, tcl = -0.2)
    }
  }
  
  
  
  if (type=="prctile"){
    bootci1=t(apply(x$mboot, 1, quantile, probs = seq(0, 1, .01),  na.rm = TRUE))
    bins1=x$BinCentres[is.na(bootci1[,1])==FALSE]
    bootci1=bootci1[is.na(bootci1[,1])==FALSE,]
    n=length(bootci1[1,])
    ## PLOT
    plot(NULL, type = "n", xlim=c(max(x$BinCentres)+w,min(x$BinCentres)-w), ylim = range(bootci1),axes=FALSE,ylab="Composite",xlab="Age",main="Percentiles")
    if(palette=="jet"){pal = colorRampPalette(rev(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")))} 
    if(palette=="BW"){
      pal = colorRampPalette(rev(c("white","black")))}
    xx=cbind(bins1,rev(bins1))
    coli=pal(50)
    for (i in 1:floor(n/2)) {
      yy <- cbind(as.vector(bootci1[,i]), rev(as.vector(bootci1[, n - i + 1])))
      polygon(xx, yy, col =coli[floor(n/2) - i + 1], border =coli[floor(n/2) - i + 1])  
    }
    for (i in c(2,11,51,91,100)){
      lines(bins1,bootci1[,i],col="grey",lty=2)
      text(min(bins1)-200,median(bootci1[1:round(length(bootci1[,1])*0.02),i]),paste(i-1,"%",sep=""),col="grey")
    }
    axis(1)
    axis(side = 1, at = seq(0, 99000, by = 500), 
         labels = FALSE, tcl = -0.2) 
    axis(2)
  }
  
  if (type=="density"){
    seqI=seq(min(na.omit(x$mboot)),max(na.omit(x$mboot)),len=1000)
    img=matrix(nrow=1000,ncol=length(x$mboot[,1]))
    
    id=seq(1,length(x$mboot[,1]),1)[rowSums(x$mboot,na.rm=T)!=0]
    id[rowSums(x$mboot,na.rm=T)!=0]
    
    for (i in seq(1,length(x$mboot[,1]),1)[rowSums(x$mboot,na.rm=T)!=0]){
      kd=density(x$mboot[i,],na.rm=T)
      img[,i]=c(approx(kd$x,kd$y,seqI)$y)
    }
    layout(matrix(c(1,1,1,2), 1, 4, byrow = TRUE))
    if(palette=="jet"){pal = colorRampPalette((c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")))} 
    if(palette=="BW"){
      pal = colorRampPalette((c("white","black")))}
    image(x$BinCentres, seqI, t(img),col = pal(100),xlab="Age",ylab="Composite",main="Density plot",axes=F, xlim=c(max(x$BinCentres)+w,min(x$BinCentres)-w))
    axis(1, cex.axis=1, xaxp=c(0,99000,99)); axis(2, cex.axis=1)
    lines(x$BinCentres,rowMeans(x$mboot, na.rm=T))
    z=matrix(1:100,nrow=1)
    x=1
    y=seq(min(img,na.rm=T),max(img,na.rm=T),len=100)
    image(x,y,z,col=pal(100),axes=FALSE,xlab="Density",ylab="")
    axis(2)
    axis(side = 1, at = seq(0, 99000, by = 500), 
         labels = FALSE, tcl = -0.2) 
  }
}
