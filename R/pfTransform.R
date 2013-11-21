pfTransform=function(IDn,
                     add=NULL,
                     Interpolate=FALSE,
                     Age=NULL,
                     method="Z-Score",
                     BasePeriod=c(-100,1e+09),
                     span=0.3,
                     RunWidth=500,
                     RunQParam=0.5,
                     stlYears=500,
                     type="BoxCox1964",
                     alpha=0.01,
                     QuantType="ALL"
){
  ## Avoid no visible binding for global variable
  paleofiresites=NULL; rm(paleofiresites)
  
  # Value for warnings
  IDChar=IDn
  
  # Check methods
  methods=c("stl", "Z-Score", "Box-Cox", "LOESS", "MinMax", "RunMed", "RunMean", "RunMin", "RunMax", "RunQuantile", "SmoothSpline", "Hurdle")
  warnmethod=method[(method %in% methods)==FALSE]
  if(length(warnmethod)!=0){stop(paste(warnmethod, "is not a valid method for pfTransform", sep=" "))}
  
  types=c("BoxCox1964", "JohnDraper")
  warntype=type[(type %in% types)==FALSE]
  if(length(warntype)!=0){stop(paste(warntype, "is not a valid type for pfBoxCox", sep=" "))}
  
  
  ## 0 Save parameters
  params=list(IDn=IDn,
              Interpolate=Interpolate,
              Age=Age,
              method=method,
              BasePeriod=BasePeriod,
              span=span,
              RunWidth=RunWidth,
              RunQParam=RunQParam,
              stlYears=stlYears,
              type=type,
              alpha=alpha)
  
  ## 1 Load charcoal paleofiredata
  if(is.null(IDn)==FALSE){
    if (is.list(IDn) & length(IDn)==2){
      
      data(paleofiredata,envir = environment())
      data(paleofiresites,envir = environment())
      
      #paleofiredata=na.omit(paleofiredata)
      # Sites are:
      IDn=IDn$SitesIDS
      # Use only paleofiredata corresponding to IDn
      paleofiredata=paleofiredata[paleofiredata[,1] %in% IDn,]
      
      ## Convert data to influx------
      if(QuantType=="INFL"){
        for(i in unique(IDn))  
          if( paleofiresites[paleofiresites$ID_SITE==i,21]!="INFL" & 
                is.na(sum(paleofiredata[paleofiredata[,1]==i,2]))==FALSE){
            temp=paleofiredata[paleofiredata[,1]==i,]
            ## Calculate Sed Acc
            d1=c()
            t1=c()
            for(k in 2:(length(temp[,1])-1)){
              d1[k]=temp[k+1,2]-temp[k-1,2]
              t1[k]=temp[k+1,3]-temp[k-1,3] 
            }
            sedacc=(d1*100)/t1
            sedacc[1]=sedacc[2]
            sedacc=c(sedacc,sedacc[length(sedacc)])
            ## Calculate Influx
            infl=(temp[,4]*sedacc)
            ## Replace in the matrix
            paleofiredata[paleofiredata[,1]==i,4]=c(infl)
          }
      }
      ##-----
      ## Add users data
      if(is.null(add)==FALSE){
        paleofiredata=rbind(paleofiredata,add$data)
        IDn=c(IDn,unique(add$data[,1]))
      }
    }
  }
  if(is.null(IDn)){
    paleofiredata=add$data
    IDn=c(unique(add$data[,1]))
  }
  
  
  
  if (is.character(IDn)){
    paleofiredata = read.csv(IDn)
    IDn=unique(paleofiredata[,1])
  }
  if (is.list(IDn) & length(IDn)>2)  {
    temp=IDn$TransData
    depths=IDn$IntDepths
    age=IDn$Age
    sites=as.numeric(colnames(temp))
    ids=matrix(nrow=length(temp[,1]),ncol=length(temp[1,]))
    for (i in 1:length(temp[,1])){
      ids[i,]=sites
    }
    ids=c(ids)
    data=c(temp)
    age=c(age) 
    depths=c(depths)
    if(length(depths)==0) depths=rep(NA,length(age))
    paleofiredata=cbind(ids,depths,age,data)
    IDn=unique( paleofiredata[,1])
  }
  if (is.matrix(IDn)){
    paleofiredata=IDn
    IDn=unique(paleofiredata[,1])
  }
  
  # 2 Interpolate TRUE
  if (Interpolate==TRUE){
    # Interpolation procedure
    if (is.null(Age)) {
      res=matrix(ncol=1,nrow=length(IDn))
      # Find the median time resolution for each paleofiredataset
      for (k in 1:length(IDn)){
        resT=diff(paleofiredata[paleofiredata[,1]==IDn[k],3])
        # Sometimes the last age is a copy of the previous one (why?)
        res[k]=c(median(resT[resT>0]))
      }
      # Find the median resolution for resampling
      step=round(median(res))
      minA=round(min(paleofiredata[,3]))
      maxA=round(max(paleofiredata[,3]))
      AgeN=seq(minA,maxA,step)
    }
    if (is.null(Ages)==FALSE) {
      AgeN=Age
      #paleofiredata=paleofiredata[paleofiredata[,3]>min(AgeN),]
      #paleofiredata=paleofiredata[paleofiredata[,3]<max(AgeN),]
      #paleofiredata=paleofiredata[paleofiredata[,1] %in% IDn,]
      IDn=unique(paleofiredata[,1])
    }
    
    # Use linear interpolation to reconstruct a matrix of raw paleofiredata
    rawI=matrix(nrow=length(AgeN),ncol=length(IDn))
    
    for (k in 1:length(IDn)){
      if(length(paleofiredata[paleofiredata[,1]==IDn[k],3])>=3){
        rawI[,k]=approx(paleofiredata[paleofiredata[,1]==IDn[k],3],paleofiredata[paleofiredata[,1]==IDn[k],4],AgeN, method = "linear")$y} else print(paste(IDChar$SiteNames[k], "has < 3 charcoal values and was excluded", sep=" "))
    }
    
    # Calculates Interpolated depths
    depthI=matrix(nrow=length(AgeN),ncol=length(IDn))
    for (k in 1:length(IDn)){
      if (is.na(sum(paleofiredata[paleofiredata[,1]==IDn[k],2]))==F){
        if(length(paleofiredata[paleofiredata[,1]==IDn[k],3])>=3){
          depthI[,k]=approx(paleofiredata[paleofiredata[,1]==IDn[k],3],paleofiredata[paleofiredata[,1]==IDn[k],2],AgeN,
                            method = "linear")$y
        }
      } else{depthI[,k]=NA}
    }
    
    ## Remove sites with less than 3 data values
    supp=c()
    for(i in 1:length(rawI[1,])){
      if(sum(!is.na(rawI[,i]))<3){supp[i]=1}else{supp[i]=0}
    }
    
    rawI=rawI[,supp==0]
    SuppSites=IDn[supp==1]
    IDn=IDn[supp==0]  
    # Space for data
    transI=matrix(nrow=length(AgeN),ncol=length(IDn))
    # Matrix of Ages (just a repeat)
    Ages=matrix(ncol=length(IDn),nrow=length(AgeN))
    for (k in 1:length(IDn)){
      Ages[,k]=c(AgeN)
    }
  }
  
  ## 3 No Interpolation:
  if (Interpolate==FALSE){
    # Which is the longest record?
    lengths=matrix(ncol=1,nrow=length(IDn))
    for (k in 1:length(IDn)){
      lengths[k]=c(length(paleofiredata[paleofiredata[,1] %in% IDn[k],1]))
    }
    m=max(lengths)
    
    # Space for paleofiredata
    transI=matrix(nrow=m,ncol=length(IDn))
    rawI=matrix(nrow=m,ncol=length(IDn))
    depthI=matrix(nrow=m,ncol=length(IDn))
    Ages=matrix(ncol=length(IDn),nrow=m)
    
    # Matrix of Ages, rawData and depths
    for (k in 1:length(IDn)){
      forNA=m-length(paleofiredata[paleofiredata[,1] %in% IDn[k],3])
      AgeTemp=c(paleofiredata[paleofiredata[,1] %in% IDn[k],3],rep(NA,forNA))
      Ages[,k]=c(AgeTemp)
      rawTemp=c(paleofiredata[paleofiredata[,1] %in% IDn[k],4],rep(NA,forNA))
      rawI[,k]=c(rawTemp)
      depthTemp=c(paleofiredata[paleofiredata[,1] %in% IDn[k],2],rep(NA,forNA))
      depthI[,k]=c(depthTemp)
    }
    ## End No Int
  }
  
  
  # Play with transformations!
  for (j in 1:length(method)){
    methodj=method[j]
    if (j>=2){rawI=transI}
    
    # Transformations
    for (k in 1:length(IDn)){
      
      tmp=cbind(Ages[,k],rawI[,k])
      tmp=na.omit(tmp)
      ## At least 3 data values! 
      if(sum(tmp[,1])>3 & IDn[k]!=882){
        # Not Tamagaucia site (882)!
        if(methodj=="stl") {
          agesI=seq(tmp[1,1],tmp[length(tmp[,1]),1],1)
          # Stl requires evenly spaced data
          forTS=approx(tmp[,1],tmp[,2],agesI)$y
          x=ts(forTS,start=1,frequency=stlYears)
          dim(x)=NULL
          stlResult=stl(x,"per")$time.series[,2]
          transI[,k]=approx(agesI,stlResult,Ages[,k])$y
        }
        if (methodj=="Z-Score") {
          mu=mean(tmp[tmp[,1]>=BasePeriod[1] & tmp[,1]<=BasePeriod[2],2])
          sigma=sd(tmp[tmp[,1]>=BasePeriod[1] & tmp[,1]<=BasePeriod[2],2])
          # No data in BasePeriod return scale:
          # Single value in BasePeriod return scale:
          if (is.na(mu) | is.na(sigma) | sigma==0){transI[,k]=approx(tmp[,1],scale(tmp[,2]),Ages[,k])$y}
          # Z-Score otherwise
          else {transI[,k]=approx(tmp[,1],(tmp[,2]-mu)/sigma,Ages[,k])$y}            
        }
        if (methodj=="Box-Cox") {
          transI[,k]=approx(tmp[,1],pfBoxCox(tmp[,2],alpha=alpha,type=type),Ages[,k])$y
        }
        if (methodj=="LOESS") {
          transI[,k]=approx(tmp[,1],predict(loess(tmp[,2]~ tmp[,1], span=span)),Ages[,k])$y
        }
        if (methodj=="MinMax") {
          transI[,k]=approx(tmp[,1],pfMinMax(tmp[,2]),Ages[,k])$y
        }
        if (methodj=="RunMed") {
          w=round(RunWidth/((max(tmp[,1])-min(tmp[,1]))/length(tmp[,1])))
          if (odd(w)) w=w else w=w+1
          transI[,k]=approx(tmp[,1],runmed(tmp[,2],w),Ages[,k])$y
        }
        if (methodj=="RunMean") {
          w=round(RunWidth/((max(tmp[,1])-min(tmp[,1]))/length(tmp[,1])))
          if (odd(w)) w=w else w=w+1
          transI[,k]=approx(tmp[,1],runmean(tmp[,2],w),Ages[,k])$y
        }
        if (methodj=="RunMin") {
          w=round(RunWidth/((max(tmp[,1])-min(tmp[,1]))/length(tmp[,1])))
          if (odd(w)) w=w else w=w+1
          transI[,k]=approx(tmp[,1],runmin(tmp[,2],w),Ages[,k])$y
        }
        if (methodj=="RunMax") {
          w=round(RunWidth/((max(tmp[,1])-min(tmp[,1]))/length(tmp[,1])))
          if (odd(w)) w=w else w=w+1
          transI[,k]=approx(tmp[,1],runmax(tmp[,2],w),Ages[,k])$y
        }
        if (methodj=="RunQuantile") {
          w=round(RunWidth/((max(tmp[,1])-min(tmp[,1]))/length(tmp[,1])))
          if (odd(w)) w=w else w=w+1
          transI[,k]=approx(tmp[,1],runquantile(tmp[,2],w,RunQParam),Ages[,k])$y
        }
        if (methodj=="SmoothSpline") {
          transI[,k]=approx(tmp[,1],smooth.spline(tmp[,1],tmp[,2],spar=span)$y,Ages[,k])$y
        }
        #         if (methodj=="GAM"){
        #           transI[,k]=approx(tmp[,1],gam(tmp[,2]~s(tmp[,1]))$fitted.values,Ages[,k])$y
        #         }
        if (methodj=="Hurdle"){
          # Transform data to count using pfMinMax
          tmp[,2]=round(pfMinMax(tmp[,2])*100)
          transI[,k]=approx(tmp[,1],hurdle(tmp[,2]~tmp[,1])$fitted.values,Ages[,k])$y
        }
      }
    }
    ## j loop end
  }
  
  ### End Return Results
  colnames(transI)=IDn
  output=structure(list(Age=structure(Ages,col.names=as.character(IDn),class="matrix" ),
                        IntDepths=structure(depthI,col.names=as.character(IDn),class="matrix" ),
                        IntData=structure(rawI,col.names=as.character(IDn),class="matrix" ),
                        TransData=structure(transI,col.names=as.character(IDn),class="matrix"),
                        Method=method,
                        params=params
  ))
  class(output)="pfTransform"  
  return(output)
  
  ###
}