pfSiteSel <- function(...) {
  
  ## Load data (bindind...)
  paleofiresites<-NULL
  data(paleofiresites,envir = environment())
  
  ## The eval function:
  theeval=function(thelist){
    eval(thelist,paleofiresites)
  }
  
  ## Retrieve all args in dots
  args=eval(substitute(alist(...)))
  
  if(length(args)>0){
    ## Eval all args
    c=lapply(args,theeval)
    d=matrix(unlist(c),ncol=length(args))
    d[is.na(d)]=FALSE
    
    ## All TRUE?
    finalTF=ifelse(rowSums(d == TRUE) == length(args), TRUE, FALSE)
    id=paleofiresites[finalTF,]$id_site
  } else id=paleofiresites$id_site
  
  ## Output:
  site_name=as.character(paleofiresites$site_name[paleofiresites$id_site %in% id])
  
  output=list(id_site=id,site_name=site_name)
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
  
  table=paleofiresites[paleofiresites$id_site %in% object$id_site,]
  rownames(table)=table$site_name
  table=subset(table, select=c("id_site", "lat", "long",
                               "elev", "min_est_age", "max_est_age", 
                               "num_dating",  "date_int", "num_samp", "l12", "rf99"))
  #print(table)
  return(table)
}

## Plot functions
plot.pfSiteSel=function(x,add=NULL,type="Map",zoom="Sites",pch="|",
                        xlim=NULL, ylim=NULL, cex=1,...)
  
{
  ## Avoid no visible binding for global variable
  paleofiresites=NULL; rm(paleofiresites)
  coast=NULL; rm(coast)
  
  data(paleofiresites,envir = environment())
  data(coast,envir = environment())
  
  ## Chronology
  if(type=="Chronology"){
    data(paleofiredata,envir = environment())
    paleofiredata=paleofiredata[paleofiredata$id_site %in% x$id_site,]
    
    IDsorted=data.frame(IDs = c(x$id_site),
                        Lat = c(paleofiresites[paleofiresites$id_site %in% x$id_site,]$lat),
                        labels = as.character(paleofiresites$site_name[paleofiresites$id_site %in% x$id_site]))
    
    IDsorted=IDsorted[with(IDsorted, order(Lat)), ]
    ## Xlim
    if(is.null(xlim)) xlim=range(paleofiredata$EST_AGE)
    
    ## Plot
    par(mar=c(4,14,2,8))
    plot(NULL, type = "n", 
         ylim = c(1,length(x$id_site)),xlim=xlim,axes=FALSE,ylab="",xlab="Age",main="Sampling resolution")
    n=c()
    for(i in 1:length(x$id_site)){
      samples=paleofiredata$EST_AGE[paleofiredata$id_site %in% IDsorted$IDs[i]]
      points(samples,rep(i,length(samples)),pch=pch,cex=cex)
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
      plot(paleofiresites$long,paleofiresites$lat,
           col="blue",xlab="Longitude",ylab="Latitude")
      lines(coast$X,coast$Y)
      points(paleofiresites[paleofiresites$id_site %in% x$id_site,]$long,
             paleofiresites[paleofiresites$id_site %in% x$id_site,]$lat, 
             bg="red",col = "red",pch = 21,xlab="Longitude",ylab="Latitude")
      if(is.null(add)==FALSE)
        points(add$metadata$LONGITUDE,
               add$metadata$LATITUDE, 
               bg="red",col = "red",pch = 21)
    }
    
    if(zoom=="Sites"|zoom=="sites"){
      # Draw map
      xl=as.vector(paleofiresites[paleofiresites$id_site %in% x$id_site,]$long)
      yl=as.vector(paleofiresites[paleofiresites$id_site %in% x$id_site,]$lat)
      
      if(is.null(xlim))
        xlim=range(xl[!is.na(xl) & is.finite(xl)])
      if(is.null(ylim))
        ylim=range(yl[!is.na(yl) & is.finite(yl)])
      
      
      plot(paleofiresites$long,paleofiresites$lat,col="blue",xlab="Longitude",ylab="Latitude",xlim=xlim,ylim=ylim)
      points(paleofiresites[paleofiresites$id_site %in% x$id_site,]$long,paleofiresites[paleofiresites$id_site %in% x$id_site,]$lat,bg="red",col = "red",pch = 21)
      lines(coast$X,coast$Y)
      if(is.null(add)==FALSE)
        points(add$metadata$LONGITUDE,
               add$metadata$LATITUDE, 
               bg="red",col = "red",pch = 21)
      
    }
  }
}

