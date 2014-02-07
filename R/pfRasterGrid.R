pfRasterGrid=function(COMP,cellsizex=NULL,cellsizey=NULL,
                      raster_extent=NULL,what=mean)
{
  # Load coast data from paleofire
  data(coast)
  
  # Create a simple table
  id_temp=pfSiteSel(as.numeric(colnames(COMP$BinnedData)))
  # Retrieve coordinates 
  sink("/dev/null");a=cbind(summary(id_temp)$LONGITUDE,
                            summary(id_temp)$LATITUDE); sink();
  
  colnames(a)=c("x","y")
  dat=data.frame(cbind(rbind(a),CHAR=c(COMP$BinnedData)))
  
  ## Then a raster
  if(is.null(raster_extent)) 
    e <- extent(c(round(range(dat$x)*2)/2,round(range(dat$y)*2)/2))
  else e <- extent(raster_extent)
  
  # Cell sizes
  if(is.null(cellsizex)) cellsizex=5 
  if(is.null(cellsizey)) cellsizey=3 
  
  # Number of rows and columns
  nc=round((e@xmax-e@xmin)/cellsizex)
  nr=round((e@ymax-e@ymin)/cellsizey)
  
  r <- raster(e, ncol=nc, nrow=nr) # Empty raster
  r1=rasterize(dat[, 1:2], r, dat[,3], fun=what) # Fill rater with mean transformed CHAR
  
  ## Basic plot
  plot(r1,col=rev(heat.colors(10)))
  lines(coast$X,coast$Y)
  points(dat[, 1:2],pch=".")
  ## Return output
  return(r1)
  
}