pfSiteSel=function(
  ID=NULL,
  Latlim=c(-360,360),
  Longlim=c(-360,360),
  Biome=NULL,
  DateInt=NULL,
  Country=NULL,
  Region=NULL,
  SiteName=NULL,
  PrefUnit=NULL,
  Elevation=c(-100000,100000),
  QuantType=NULL,
  L12=NULL,
  RF99=NULL
)
{
  #Parameters and data
  data(paleofiresites,envir = environment())
  coln=length(paleofiresites[1,])
  
  ## ID_SITE
  if(is.null(ID)){ID=unique(paleofiresites$ID_SITE)}
  ## SiteName
  if(is.null(SiteName)){SiteName=unique(paleofiresites$SITE_NAME)}
  ## Country
  if(is.null(Country)){Country=unique(paleofiresites$ID_COUNTRY)}
  ## PrefUnit
  if(is.null(PrefUnit)){PrefUnit=unique(paleofiresites$PREF_UNITS)}
  ## PrefUnit
  if(is.null(Region)){Region=unique(paleofiresites$ID_REGION)}
  ## Biome
  if(is.null(Biome)){Biome=unique(paleofiresites$BIOME)}
  ## QTYPE
  if(is.null(QuantType)){QuantType=unique(paleofiresites$QTYPE)}
  ## RF99
  if(is.null(RF99)){RF99=unique(paleofiresites$RF99)}
  ## L12
  if(is.null(L12)){L12=unique(paleofiresites$L12)}
  # Date Numbers
  paleofiresites=cbind(paleofiresites,(paleofiresites$MAX_EST_AGE-paleofiresites$MIN_EST_AGE)/paleofiresites$NUM_DATING)
  colnames(paleofiresites)[coln+1] <- "DATE_MEAN"
  
  paleofiresites[!is.finite(paleofiresites[,coln+1]),coln+1]=NA
  if(is.null(DateInt)){DateInt=NA}
  
  
  # Selection process
  
  IDn=paleofiresites[  paleofiresites$ID_SITE %in% ID 
                       & paleofiresites$LONGITUDE>=Longlim[1] 
                       & paleofiresites$LONGITUDE<=Longlim[2] 
                       & paleofiresites$LATITUDE>=Latlim[1] 
                       & paleofiresites$LATITUDE<=Latlim[2]  
                       & paleofiresites$ELEV>=Elevation[1]
                       & paleofiresites$ELEV<=Elevation[2]
                       & paleofiresites$BIOME %in%  Biome 
                       & paleofiresites$ID_COUNTRY %in%  Country 
                       & paleofiresites$RF99 %in%  RF99
                       & paleofiresites$L12 %in%  L12
                       & paleofiresites$ID_REGION %in% Region
                       & paleofiresites$PREF_UNITS %in% PrefUnit
                       & paleofiresites$QTYPE %in% QuantType 
                       & paleofiresites$SITE_NAME %in% SiteName, 1]
  
  if(is.na(DateInt)==FALSE)
  IDn=paleofiresites[paleofiresites$ID_SITE %in% IDn & 
                   paleofiresites$DATE_MEAN<DateInt,1] 
  
  #& paleofiresites$ID_COUNTRY %in% Country 
  
  # Returns NA (Why? to be checked)
  
  IDn=as.numeric(na.omit(IDn))
  
  # Site Names
  SiteNames=as.character(paleofiresites$SITE_NAME[paleofiresites$ID_SITE %in% IDn])
  #
  output=list(SitesIDS=IDn,SiteNames=SiteNames)
  class(output)="pfSiteSel"
  return(output)
  
}

## Summary function
summary.pfSiteSel=function(object,...){
  ## Avoid no visible binding for global variable
  paleofiresites=NULL; rm(paleofiresites)
  paleofiredata=NULL; rm(paleofiredata)
  
  data(paleofiresites,envir = environment())
  data(paleofiredata,envir = environment())
  coln=length(paleofiresites[1,])
  
  NUM_SAMP=c()
  for (i in 1:length(object$SitesIDS)){
    NUM_SAMP[i]=c(length(paleofiredata[paleofiredata[,1]==object$SitesIDS[i],1]))
  }
  table=cbind(paleofiresites[paleofiresites$ID_SITE %in% object$SitesIDS,],NUM_SAMP)
  rownames(table)=table$SITE_NAME
  table=subset(table, select=c(1,3,4,5,16,17,18,coln+1))
  print(table)
  return(table)
}

## Plot functions
plot.pfSiteSel=function(x,type="Map",zoom="Sites",...){
  
  ## Avoid no visible binding for global variable
  paleofiresites=NULL; rm(paleofiresites)
  coast=NULL; rm(coast)
  
  data(paleofiresites,envir = environment())
  data(coast,envir = environment())
  
  ## Chronology
  if(type=="Chronology"){
    data(paleofiredata,envir = environment())
    paleofiredata=paleofiredata[paleofiredata$ID_SITE %in% x$SitesIDS,]
    
    IDsorted=data.frame(IDs = c(x$SitesIDS),
                        Lat = c(paleofiresites[paleofiresites$ID_SITE %in% x$SitesIDS,]$LATITUDE),
                        labels = as.character(paleofiresites$SITE_NAME[paleofiresites$ID_SITE %in% x$SitesIDS]))
    
    IDsorted=IDsorted[with(IDsorted, order(Lat)), ]
    
    par(mar=c(4,14,2,8))
    plot(NULL, type = "n", xlim=c(max(paleofiredata$EST_AGE),min(paleofiredata$EST_AGE)), ylim = c(1,length(x$SitesIDS)),axes=FALSE,ylab="",xlab="Age",main="Sampling resolution")
    n=c()
    for(i in 1:length(x$SitesIDS)){
      samples=paleofiredata$EST_AGE[paleofiredata$ID_SITE %in% IDsorted$IDs[i]]
      points(samples,rep(i,length(samples)),pch="|")
      n[i]=length(samples)
    }
    axis(2, at=seq(1,length(IDsorted$IDs),1), labels = FALSE)   
    IDsorted$labels=gsub("[\x87]", "c", IDsorted$labels) 
    IDsorted$labels=gsub("[\x85]", "a", IDsorted$labels) 
    IDsorted$labels=gsub("[\x82]", "e", IDsorted$labels) 
    IDsorted$labels=gsub("[\x8a]", "e", IDsorted$labels) 
    text(y = seq(1,length(IDsorted$IDs),1), par("usr")[1], labels = IDsorted$labels, srt = 0, pos = 2, xpd = TRUE)
    axis(side = 1, at = seq(0, 99000, by = 500), labels = FALSE, tcl = -0.2) 
    axis(4, at=seq(1,length(IDsorted$IDs),1), labels = FALSE)    
    text(y = seq(1,length(IDsorted$IDs),1), par("usr")[2], labels = paste(round(IDsorted$Lat,digits=1),"/",n,sep=""), srt = 0, pos = 4, xpd = TRUE)
    
    paste(round(IDsorted$Lat,digits=1),"/",n,sep="")
    axis(1)
  }
  
  ## MAPS
  if(type=="Map"){
    
    if(zoom=="World"|zoom=="world"){
      plot(paleofiresites$LONGITUDE,paleofiresites$LATITUDE,col="blue",xlab="Longitude",ylab="Latitude")
      lines(coast$X,coast$Y)
      points(paleofiresites[paleofiresites$ID_SITE %in% x$SitesIDS,]$LONGITUDE,paleofiresites[paleofiresites$ID_SITE %in% x$SitesIDS,]$LATITUDE, bg="red",col = "red",pch = 21,xlab="Longitude",ylab="Latitude")
    }
    
    if(zoom=="Sites"|zoom=="sites"){
      # Draw map
      xl=as.vector(paleofiresites[paleofiresites$ID_SITE %in% x$SitesIDS,]$LONGITUDE)
      yl=as.vector(paleofiresites[paleofiresites$ID_SITE %in% x$SitesIDS,]$LATITUDE)
      
      xlim=range(xl[!is.na(xl) & is.finite(xl)])
      ylim=range(yl[!is.na(yl) & is.finite(yl)])
    
      plot(paleofiresites$LONGITUDE,paleofiresites$LATITUDE,col="blue",xlab="Longitude",ylab="Latitude",xlim=xlim,ylim=ylim)
      points(paleofiresites[paleofiresites$ID_SITE %in% x$SitesIDS,]$LONGITUDE,paleofiresites[paleofiresites$ID_SITE %in% x$SitesIDS,]$LATITUDE,bg="red",col = "red",pch = 21)
      lines(coast$X,coast$Y)}
  }
}

