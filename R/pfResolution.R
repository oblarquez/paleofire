pfResolution=function(ID,AgeLim=NULL){
  
  # Temporal resolution of records
  
  paleofiredata=NULL; rm(paleofiredata)
  
  data(paleofiredata, envir = environment())
  IDn=ID$SitesIDS
  # Use only paleofiredata corresponding to IDn
  if (is.null(AgeLim)){
  paleofiredata=paleofiredata[paleofiredata[,1] %in% IDn,]
  } else {
    paleofiredata=paleofiredata[paleofiredata[,1] %in% IDn,]
    paleofiredata=paleofiredata[paleofiredata$EST_AGE>AgeLim[1] & paleofiredata$EST_AGE<AgeLim[2],]
  }
  
  meanres=c()
  medianres=c()
  sdres=c()
  
  for (i in 1:length(IDn)){
    meanres[i]=mean(diff(paleofiredata[paleofiredata$ID_SITE %in% IDn[i],]$EST_AGE))
    medianres[i]=median(diff(paleofiredata[paleofiredata$ID_SITE %in% IDn[i],]$EST_AGE))
    sdres[i]=sd(diff(paleofiredata[paleofiredata$ID_SITE %in% IDn[i],]$EST_AGE))
  }
  
  res=data.frame(ID_SITE=as.numeric(ID$SitesIDS),
                 SITE_NAME=ID$SiteNames,
                 MeanRes=as.numeric(meanres),
                 MedianRes=as.numeric(medianres),
                 SdRes=as.numeric(sdres))
  return(res)
  
}