pfSiteSelM <- function(...) {
  paleofiresites<-NULL
  data(paleofiresites,envir = environment())
  
  theeval=function(thelist){
    eval(thelist,paleofiresites)
  }
  
  args=eval(substitute(alist(...)))
  if(length(args)>0){
  c=lapply(args,theeval)
  d=matrix(unlist(c),ncol=length(args))
  finalTF=ifelse(rowSums(d == TRUE) == length(args), TRUE, FALSE)
  id=paleofiresites[finalTF,]$ID_SITE
  } else id=paleofiresites$ID_SITE
  
  SiteNames=as.character(paleofiresites$SITE_NAME[paleofiresites$ID_SITE %in% id])
  #
  output=list(SitesIDS=id,SitesNames=SiteNames)
  class(output)="pfSiteSel"
  return(output)
}


