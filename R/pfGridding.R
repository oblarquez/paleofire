## 
pfGridding=function(data,cell_sizex=NULL,
                    cell_sizey=NULL,
                    age=0,
                    cell_size=NULL,
                    time_buffer=NULL,
                    distance_buffer=NULL,
                    threshold=0.5,
                    raster_extent=NULL,
                    elevation_range=NULL,
                    proj4=NULL,
                    sea_mask=FALSE,
                    verbose=TRUE){
  
  ## pfTransform object
  if(class(data)=="pfTransform"){
    sink("/dev/null")
      data<-data.frame(
        x=rep(summary(data$params$ID)$LONGITUDE,each=length(data$TransData[,1])),
        y=rep(summary(data$params$ID)$LATITUDE,each=length(data$TransData[,1])),
        age=c(data$Age),
        char=c(data$TransData));
    sink();
    data=na.omit(data)
    ## Remove extreme outliers i.e. 3*sd (sometimes happens... why? probably baseperiod related)
    data=data[data[,4]<3*sd(data[,4]) & data[,4]>-(3*sd(data[,4])) ,]
  }
  
  
  ## Limit the input to age+-time_buffer
  if(is.null(time_buffer)) time_buffer=500
  data=data[data[,3]>=age-time_buffer & data[,3]<=age+time_buffer,]
  
  
  # source("/Users/Olivier/Documents/BorealTreeCover/final/triCube.R")
  ## Load countries with lakes from http://www.naturalearthdata.com/downloads/10m-cultural-vectors/
  # load(file="/Users/Olivier/Documents/BorealTreeCover/final/world_map.rda")
  
  
  if(is.null(proj4))
    proj4<-"+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  
  xy <-cbind(data[,1],data[,2])
  
  ### LOAD DEM if required
  if(is.null(elevation_range)==FALSE | sea_mask==TRUE){
    cat("Preparing data and loading DEM...")
    cat("\n")
    ## dem is GMTED2010
    data(dem)
    if(is.null(raster_extent)) v=extent(xy) else v=extent(raster_extent)
    dem=crop(dem,v)
    dem1 <- projectRaster(dem, crs=proj4)
    #plot(dem1)
  }
  ###
  
  
  dat2=data.frame(project(xy, proj4))
  colnames(dat2)=c("x","y")
  if(is.null(raster_extent)) {
    e <- extent(c(round(range(dat2$x)*2)/2,round(range(dat2$y)*2)/2))
  } else {
    toto=cbind(x=c(raster_extent[c(1,1,2,2)], raster_extent[1]+(raster_extent[2]-raster_extent[1])/2),
               y=c(raster_extent[c(3,4,3,4)], raster_extent[3]))
    e <- extent(project(toto, proj4))
  }
  
  
  # Cell sizes
  if(is.null(cell_size)){
    if(is.null(cell_sizex)) cell_sizex=200000
    if(is.null(cell_sizey)) cell_sizey=200000
  } else {
    cell_sizex=cell_size
    cell_sizey=cell_size
  }
  # Number of rows and columns
  nc=ceiling((e@xmax-e@xmin)/cell_sizex)
  nr=ceiling((e@ymax-e@ymin)/cell_sizey)
  
  r <- raster(e, ncol=nc, nrow=nr) # Empty raster
  projection(r)<-proj4
  
  if(is.null(distance_buffer)) distance_buffer=300000
  
  ## Elevation stuff (median elevation in each predicted cell)  
  if(is.null(elevation_range)==FALSE | sea_mask==TRUE){
    temp=rasterToPoints(dem1)
    temp1=rasterize(temp[, 1:2], r, temp[,3], fun=mean)
    z=raster::intersect(temp1,r)
    # plot(temp1)
    elev1=rasterToPoints(temp1)[,3]
    dat1=as.data.frame(rasterToPoints(z))
  } else dat1=as.data.frame(rasterToPoints(r))
  
  
  ## Initiate a percentage counter
  if(verbose==TRUE){
    percent=seq(10,100,by=10)
    values=round(percent*length(dat1[,1])/100)
    cat("Spatio-temporal weighted interpolation...")
    cat("\n")
    cat("Percentage done: ")
  }
  
  #ptm <- proc.time()
  
  ## Main loop time--distance weighting 
  dat1[,3]=c()
  for(i in 1:length(dat1[,1])){
    d=pointDistance(dat1[i,1:2], dat2,longlat=FALSE)
    d=cbind(dat2,data[,3],d,triCube(d,distance_buffer))
    
    ## Elevation range search
    if(is.null(elevation_range)==FALSE ){
      elev2=extract(dem1,d[,1:2])
      elev1=extract(temp1,dat1[i,1:2])
      elev2[elev2<elev1-elevation_range]=NA
      elev2[elev2>elev1+elevation_range]=NA
      d=cbind(d,elev2)
      d[is.na(d[,6]),5]=0
    }
    d1=d[d[,5]!=0,]
    
    ## Time weight
    d1=cbind(d1,triCube(age-d1[,3],time_buffer))
    # head(d)
    if(is.null(elevation_range)==FALSE){
      colnames(d1)=c("x","y","age","dist","dweight","elev","tweight")
    } else colnames(d1)=c("x","y","age","dist","dweight","tweight")
    
    d1$weight=d1$dweight*d1$tweight
    d1$data=data[d[,5]!=0,4]
    
    ## Combine time and distance weights
    if(sum(d1$weight,na.rm=TRUE)>threshold){
      dat1[i,3]=sum(d1$data*d1$weight,na.rm=TRUE)/sum(d1$weight,na.rm=TRUE)
    } else dat1[i,3]=NA
    
    ## Verbose
    if(i %in% values & verbose==TRUE)
      cat(percent[values==i]," ",sep="")
    #cat(i," ")
  }
  
  r1=rasterize(dat1[, 1:2], r, dat1[,3], fun=mean) # Fill raster with mean value 
  #   plot(r1)
  
  
  ## SEA MASK  
  if(sea_mask==TRUE){
    r2 <- is.na(temp1)  ## sea is now 1
    r2[r2==1]=NA
    #plot(r2) 
    # !!!!!!MASK
    r1=(r1-r2)
  }
  
  # plot(r1-r2)
  
  dat1=as.data.frame(rasterToPoints(r1))
  
  ## End of DATA part --------------------------------------------------------------------
  
  out=list(raster=r1,df=dat1,proj4=proj4,extent=e,points=dat2)
  class(out)="pfGridding"
  return(out)
  
}

plot.pfGridding=function(x,continuous=TRUE,
                col_class=NULL,
                col_lim=NULL,
                xlim=NULL,ylim=NULL,empty_space=10,
                cpal="YlGn",
                anomalies=TRUE,
                file=NULL,points=FALSE,...){
  
  x$df=as.data.frame(rasterToPoints(x$raster))
  # Define classes for colors
  if(is.null(col_class)){
    if(anomalies==TRUE){
      b1=seq(floor(min(x$df$layer,na.rm=TRUE)),0,len=5)
      b2=seq(0,ceiling(max(x$df$layer,na.rm=TRUE)),len=5)
      breaks=c(b1,b2[2:5])
      breaks=unique(breaks)
    } else  breaks=seq(floor(min(x$df$layer,na.rm=TRUE)),ceiling(max(x$df$layer,na.rm=TRUE)),len=10)
  }
  if(is.numeric(col_class) & length(col_class)==1){
    if(anomalies==TRUE){
      b2=seq(0,ceiling(max(abs(x$df$layer,na.rm=TRUE)))+col_class,by=col_class)
      breaks=c(-rev(b2),b2[2:length(b2)])
    } else breaks=seq(floor(min(x$df$layer,na.rm=TRUE)),ceiling(max(x$df$layer,na.rm=TRUE))+col_class,by=col_class)
  }
  if(is.numeric(col_class) & length(col_class)>1){
    breaks=col_class
  }
  # End of options
  
  ## Define color limils
  if(is.null(col_lim))
    col_lim=c(floor(min(x$df$layer,na.rm=TRUE)),ceiling(max(x$df$layer,na.rm=TRUE)))
  
  x$df=cbind(x$df,class=cut(x$df$layer,breaks))
  #   data=na.omit(data)
  
  data(coast)
  coast=coast[coast$X<=180,]
  # plot((coast$X),(coast$Y),type="l")  
  xy=cbind(coast[,2],coast[,1])
  coast=data.frame(project(xy, x$proj4))
  colnames(coast)=c("x","y")
  # plot(round(coast$x),round(coast$y),type="l")
  
  ## LIMITS
  if(is.null(xlim)){
    xplus=(x$extent@xmax-x$extent@xmin)*empty_space/100
    xlim=c(x$extent@xmin-xplus,x$extent@xmax+xplus)}
  
  if(is.null(ylim)){
    yplus=(x$extent@ymax-x$extent@ymin)*empty_space/100
    ylim=c(x$extent@ymin-yplus,x$extent@ymax+yplus)}
  
  ## Crop coast using limits
    coast=coast[coast$x>xlim[1]-8000000 & coast$x<xlim[2]+8000000 &
                  coast$y>ylim[1]-8000000 & coast$y<ylim[2]+8000000,]
  #plot(coast[,1],coast[,2],type="l")
  
  x$points=data.frame(na.omit(x$points))
  colnames(x$points)=c("x","y")
  
  if(anomalies==TRUE) {
    cpal="RdBu"
    pal=rev(brewer.pal(9,cpal))} else pal=brewer.pal(9,cpal)
  
  testcol  <- colorRampPalette(pal)
  
  ## LIMITS
  if(is.null(xlim)){
    xplus=(x$extent@xmax-x$extent@xmin)*0.1
    xlim=c(x$extent@xmin-xplus,x$extent@xmax+xplus)}
  
  if(is.null(ylim)){
    yplus=(x$extent@ymax-x$extent@ymin)*0.1
    ylim=c(x$extent@ymin-yplus,x$extent@ymax+yplus)}
  
  ## SAME COLORS 
  if(continuous==FALSE & is.numeric(col_class) & length(col_class)>1){
    c1=breaks+(mean(diff(breaks)/2))
    c3=cbind(xlim[2]+1e+6,ylim[2]+1e+6,c1)
    c3=c3[1:(length(c3[,1])-1),]
    c3=data.frame(c3,rr=(cut(c3[,3],breaks)))
    colnames(c3)=colnames(x$df)
    x$df=rbind(x$df,c3)
  }
  ## 
  pale=testcol(length(unique(x$df$class)))
  
  #display.brewer.pal(12,"Spectral")
  #pal=c(pal[9],pal[8],pal[6:1])
  ## On fait une carte avec ggplot2
  #pdf(file="/Users/Olivier/Desktop/x$dfMap.pdf",height=6,width=9)
  if(continuous==FALSE){
    p=ggplot(x$df) +
      geom_polygon(data=coast,aes(x=x,y=y),colour="grey80",fill="grey80")+
      geom_raster(data=x$df,aes(x, y, fill = class))+
      scale_fill_manual(values = pale,name="")+
      coord_cartesian(xlim=xlim,ylim=ylim)+xlab("Longitude")+ylab("Latitude")+
      theme_bw(base_size = 16)
    if(points==TRUE) p=p+geom_point(data=x$points,aes(x,y),colour="grey40")

  } else {
    p=ggplot(x$df) +
      geom_polygon(data=coast,aes(x=x,y=y),colour="grey80",fill="grey80")+
      geom_raster(data=x$df,aes(x, y, fill = layer))+
      scale_fill_gradient2(high=pal[9],low=pal[1],mid="white",limits=col_lim)+
      coord_cartesian(xlim=xlim,ylim=ylim)+xlab("Longitude")+ylab("Latitude")+
      theme_bw(base_size = 16)
    if(points==TRUE) p=p+geom_point(data=x$points,aes(x,y),colour="grey40")
      
  }
  p
  if(is.null(file)==FALSE){
    projection(x$raster)<-CRS(x$proj4)
    writeRaster(x$raster, filename=file, format="GTiff",overwrite=TRUE)
  }
  return(p)
}

## -------------------------------------------------------------------------------------------
#'   Tukey's Tricube weight function
#'   
#'   From the EGRET package http://usgs-r.github.io/EGRET/
#'   Robert Hirsch and Laura De Cicco
#'   
#'      Computes the tricube weight function on a vector of distance values (d),
#'      based on a half-window width of h,
#'      and returns a vector of weights that range from zero to 1.
#'
#' @param d numeric vector of distances from the point of estimation to the given sample value
#' @param h numeric value, the half-window width, measured in the same units as d
#' @keywords statistics weighting
#' @return w numeric vector of weights, all 0<=w<=1
#' @export
#' @examples
#'  h<-10
#'  d<-c(-11,-10,-5,-1,-0.01,0,5,9.9,10,20)
#'  triCube(d,h)
triCube<-function(d,h) {
  #  triCube is Tukey tricubed weight function
  #    first argument, d, is a vector of the distances between the observations and the estimation point
  #    second argument, h, is the half window width
  #    it returns a vector of weights (w) for the observations in the vector, d
  n<-length(d)
  zero<-rep(0,n)
  ad<-abs(d)
  w<-(1-(ad/h)^3)^3
  w<-pmax(zero,w)
  return(w)
}