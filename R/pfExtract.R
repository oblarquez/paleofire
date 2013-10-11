pfExtract=function(IDn){
  ## Avoid no visible binding for global variable
  paleofiredata=NULL; rm(paleofiredata)
  
  # Extract data for sites
  data(paleofiredata)
  if(is.numeric(IDn)){  Ext=paleofiredata[paleofiredata[,1] %in% IDn,] } else Ext=paleofiredata[paleofiredata[,1] %in% IDn$SitesIDS,]
  rm(paleofiredata,envir = globalenv())
  return(Ext)
}