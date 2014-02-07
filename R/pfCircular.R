pfCircular=function(comp,b=NULL,conf=c(0.05,0.95),nboot=1000,AgeLim=NULL){
  
  ## R function developped from SEA.m   
  
  
  ## Load matrix
  T=comp$BinnedData
  
  ## Define Age limits
  if(is.null(AgeLim)==FALSE){
    T=T[comp$BinCentres>AgeLim[1] & comp$BinCentres<AgeLim[2], ]
  }
  
  ## Block size calculus
  if(is.null(b)==TRUE){
    b=c()
    for(i in 1:length(T[1,])){
      r=cor(T[1:length(T[,1])-1,i],T[2:length(T[,1]),i],use="pairwise.complete.obs")
      yb=2*(-1/log(abs(r)))
      b[i]=c(ceiling(yb/mean(diff(comp$BinCentres))))
    }
    b[b==0 | b==1 | is.na(b) | is.finite(b)==FALSE]=2} else {b=rep(b,length(T[1,]))}
  
  
  ## Arrange data
  a=matrix(nrow=b,ncol=length(T[1,]))
  a[is.na(a)]=-999
  T_=rbind(a,T,a)
  
  ## Declare values for the boot process
  y_m=matrix(ncol=nboot,nrow=length(T[,1]))
  cat("# of Bootstrap:")
  
  for(k in 1:nboot){
    y_n=matrix(nrow=length(T[,1]),ncol=length(T[1,]))
    for(i in 1:length(T[1,])){
      n=ceiling(length(T_[,1])/b[i])
      q=trunc(length(T_[,1])-b[i])+1
      y=matrix(nrow=b[i],ncol=q)
      for(j in 1:q){
        y[,j]=c(T_[seq((j-1)+1,(j-1)+b[i],1),i])
      }
      o=sample(1:q,n*2,replace=TRUE)
      yy=c(y[,o])
      yy=yy[!yy == -999] 
      y_n[,i]=yy[1:length(T[,1])]
    }
    y_m[,k]=rowMeans(y_n,na.rm=TRUE)
    
    if(k %in% seq(0,nboot,10)) cat("", k)
  }
  
  ## Compile conf intervals
  boots=t(apply(y_m, 1, quantile, probs = conf,  na.rm = TRUE))
  
  ## Values for output
  if(is.null(AgeLim)==FALSE){
    yr=comp$BinCentres[comp$BinCentres>AgeLim[1] & comp$BinCentres<AgeLim[2]]
    Ci=comp$BootCi[comp$BinCentres>AgeLim[1] & comp$BinCentres<AgeLim[2], ]
    BootMean=comp$BootMean[comp$BinCentres>AgeLim[1] & comp$BinCentres<AgeLim[2]]
  } else {
    yr=comp$BinCentres
    Ci=comp$BootCi
    BootMean=comp$BootMean
  }
  
  ## Output
  output=structure(list(BootCirc=structure(boots,row.names = as.character(yr),col.names=conf, class = "matrix"),
                        conf=conf,
                        yr=yr,
                        BootCi=Ci,
                        BootMean=BootMean))
  class(output)="pfCircular"
  return(output)  
}

##--------------------------------------------------------------------------------------------------------
plot.pfCircular=function(x,...){
  ## Plot
  
  t=c(x$BootMean,x$BootCirc)
  plot(x$yr,x$BootMean,type="o",
       ylim=c(min(t,na.rm=TRUE),max(t,na.rm=TRUE)),
       xlim=c(max(x$yr),min(x$yr)),xlab="Age (cal yr BP)",ylab="Composite",font.main=1, lab=c(8,5,5), 
       cex.lab=1, pch=16, cex=0.5,axes=F, mgp=c(2,0,0))
  for (i in 1:length(x$BootCirc[1,])){
    lines(x$yr,x$BootCirc[,i],lty=2)
    text(min(x$yr)-200,x$BootCirc[1,i],paste(x$conf[i]*100,"%",sep=""),col="black")
  }
  axis(1); axis(2, cex.axis=1)
  axis(side = 1, at = seq(0, 99000, by = 500), 
       labels = FALSE, tcl = -0.2)  
}




