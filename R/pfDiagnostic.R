pfDiagnostic=function(ID,
                      add=NULL,
                      Age=0,
                      Interpolate=FALSE,
                      method="Box-Cox",
                      BasePeriod=c(-100,1e+09),
                      span=0.3,
                      RunWidth=500,
                      RunQParam=0.5,
                      stlYears=500,
                      alpha=0.01,
                      type="BoxCox1964",
                      FileName="Diagnostic.pdf",
                      QuantType="ALL"
){
  
  ## Avoid no visible binding for global variable
  coast=NULL; rm(coast)
  
  ## Load data
  if (is.list(ID) & length(ID)==2){
    data(paleofiredata,envir = environment())
    data(paleofiresites,envir = environment())
    #paleofiredata=na.omit(paleofiredata)
    ID=ID$id_site
    # Use only paleofiredata corresponding to ID
    paleofiredata=paleofiredata[paleofiredata[,1] %in% ID,]
    
    ## Convert data to influx------
#     if(QuantType=="INFL"){
#       for(i in unique(ID))  
#         if( paleofiresites[paleofiresites$id_site==i,]$QTYPE!="INFL" & 
#               is.na(sum(paleofiredata[paleofiredata[,1]==i,2]))==FALSE){
#           temp=paleofiredata[paleofiredata[,1]==i,]
#           ## Calculate Sed Acc
#           d1=c()
#           t1=c()
#           for(k in 2:(length(temp[,1])-1)){
#             d1[k]=temp[k+1,2]-temp[k-1,2]
#             t1[k]=temp[k+1,3]-temp[k-1,3] 
#           }
#           sedacc=(d1*100)/t1
#           sedacc[1]=sedacc[2]
#           sedacc=c(sedacc,sedacc[length(sedacc)])
#           ## Calculate Influx
#           infl=(temp[,4]*sedacc)
#           ## Replace in the matrix
#           paleofiredata[paleofiredata[,1]==i,4]=c(infl)
#         }
#     }
    ##-----
    ## Add users data
    if(is.null(add)==FALSE){
      paleofiredata=rbind(paleofiredata,add$data)
      ID=c(ID,unique(add$data[,1]))
      paleofiresites=merge(paleofiresites,add$metadata,all=TRUE)
    }
  }
  
  if (is.character(ID)){
    paleofiredata = read.csv(ID)
    ID=unique( paleofiredata[,1])
  }
  
  # Warnings for the age depths model
  warn=matrix(nrow=length(ID),ncol=1)
  
  
  
  # plot the original data, histograms, and transformed data -- use one of the following and see dev.off() below
  pdf(file=FileName,width=8.5, height=11.0,onefile=T)  # create a .pdf file of plots
  
  layout(matrix(c(1,1,2,2,3,4,5,6),4,2,byrow=TRUE))
  
  for (k in 1:length(ID)){
    # Sites identifier
    ID1=ID[k]
    # Site info
    siteinfo=paleofiresites[paleofiresites[,1]==ID1,]
    # List for pfTransform function
    IDt=list(id_site=ID1,site_name=as.character(siteinfo$site_name))
    
    
    # Little TRICK for sites without depths = Assume depths == 1:length(Record)
    # i.e. continuous sampling
    if (is.na(sum(paleofiredata[paleofiredata[,1]==ID1,2]))==TRUE){  
      depth=c(1:length(paleofiredata[paleofiredata[,1]==ID1,2]))
      warn[k,]=1
    }else{warn[k,]=0
          depth=paleofiredata[paleofiredata[,1]==ID1,2]}
    
    
    # plot raw data
    #par(fig=c(0,1,.70,.95), mar=c(3.1, 4.1, 2.1, 4.1))
    plot(paleofiredata[paleofiredata[,1]==ID1,3],paleofiredata[paleofiredata[,1]==ID1,4], xlim=c(20000,-100), axes=F, mgp=c(2,0,0),
         font.main=1, lab=c(10,5,5), 
         ylab=paste("Char", siteinfo$pref_units,sep=" "),xlab=" " ,cex.lab=0.8, pch=16, cex=0.5, type="o")
    axis(1, cex.axis=0.8, xaxp=c(0,22000,22)); axis(2, cex.axis=0.8); axis(4, cex.axis=0.8)
    title(paste(siteinfo$site_name, ", Site ID: ", siteinfo$id_site, sep=""),cex.main=1.5)
    
    # Transform data
    if(ID1<10000) {tr=IDt 
                   add1=NULL } 
    if(ID1>10000){tr=NULL 
                  add1=c()
                  add1$data=paleofiredata[paleofiredata$id_site %in% ID1,]}
    

      tr=pfTransform(ID=tr,add=add1,
                     Interpolate=Interpolate,
                     method=method,
                     stlYears=stlYears,
                     Age=Age,
                     BasePeriod=BasePeriod,
                     span=span,
                     RunWidth=RunWidth,
                     RunQParam=RunQParam,
                     type=type,
                     alpha=alpha)
    
    # plot transformed data
    if(sum(is.na(tr$TransData))!=length(tr$TransData)){
      #par(new=T, fig=c(0,1,.4,.65), mar=c(3.1, 4.1, 2.1, 4.1))
      plot(tr$Age, tr$TransData,xlim=c(20000,-100), axes=F, mgp=c(2,0,0), lab=c(10,5,5),
           ylab=paste(method,"transformed data",sep=" "), xlab="Age", cex.lab=0.8, pch=16, cex=0.5, type="o")
    axis(1, cex.axis=0.8, xaxp=c(0,22000,22)); axis(2, cex.axis=0.8); axis(4, cex.axis=0.8)
    } else plot.new()
    
    # plot histograms of raw and transformed data, and likelhood profile
    #par(new=T, fig=c(0,.4,.0,.35), mar=c(7.1, 4.1, 2.1, 2.1))
    hist(paleofiredata[paleofiredata[,1]==ID1,4], xlab=paste("Char", siteinfo$pref_units,sep=" "), cex.lab=0.8, cex.axis=0.8, mgp=c(2,1,0), main=NULL)
    
    if(sum(is.na(tr$TransData))!=length(tr$TransData)){
      #par(new=T, fig=c(.3,.7,.0,.35), mar=c(7.1, 4.1, 1.1, 2.1))
      hist(tr$TransData, xlab=paste(method,"transformed data",sep=" "), cex.lab=0.8, cex.axis=0.8, mgp=c(2,1,0), main=NULL)
    } else plot.new()
    # MAP
    #data(paleofiresites)
    data(coast, envir = environment())  
    # Draw map
    plot(paleofiresites$long,paleofiresites$lat,col="blue",
         xlim=c(paleofiresites[paleofiresites$id_site %in% IDt,]$long-20,paleofiresites[paleofiresites$id_site %in% IDt,]$long+20),  
         ylim=c(paleofiresites[paleofiresites$id_site %in% IDt,]$lat-20,paleofiresites[paleofiresites$id_site %in% IDt,]$lat+20),
         xlab="Longitude",ylab="Latitude")
    points(paleofiresites[paleofiresites$id_site %in% IDt,]$long,paleofiresites[paleofiresites$id_site %in% IDt,]$lat,
           bg="red",col = "red",pch = 21)
    lines(coast$X,coast$Y)
    
    # Age depth
    plot(paleofiredata[paleofiredata[,1]==ID1,3],depth,type="l",
         xlab="Age (yrs BP)",ylab="Depth (m)",xlim=c(max(paleofiredata[paleofiredata[,1]==ID1,3]),min(paleofiredata[paleofiredata[,1]==ID1,3])),
         ylim=c(max(depth),min(depth)),
         , cex.lab=0.8, cex.axis=0.8, mgp=c(2,1,0), main=NULL);
    axis(1, cex.axis=0.8); axis(2, cex.axis=0.8); axis(4, cex.axis=0.8)
    text(y=min(depth), x=mean(paleofiredata[paleofiredata[,1]==ID1,3]), labels = paste("# Dates=", siteinfo$num_dating,sep=""))
    if (warn[k,]==1){
      text(y=mean(depth), x=mean(paleofiredata[paleofiredata[,1]==ID1,3]), labels = "Warning depths assumed continuous",col="red")  
    }
    
  }
  
  dev.off()
}


