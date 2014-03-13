pfExtract=function(ID){
  ## Avoid no visible binding for global variable
  paleofiredata=NULL; rm(paleofiredata)
  
  # Extract data for sites
  data(paleofiredata, envir = environment())
  if(is.numeric(ID)){  Ext=paleofiredata[paleofiredata[,1] %in% ID,] } else Ext=paleofiredata[paleofiredata[,1] %in% ID$id_site,]
  return(Ext)
}