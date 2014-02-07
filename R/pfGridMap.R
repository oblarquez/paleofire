pfGridMap = function(PF, bins, id=NULL, site.dat=NULL,
                     fig.base.name=NULL, base.map='coasts',
                     grd.res=NULL, grd.ext=c(-180,180,-90,90), 
                     proj.map=TRUE, proj.str="+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs", 
                     max.kr.dist=NULL, allow.interpolate=FALSE
) {
  
  data(countriesCoarse)  # A dataset in rworldmap used in the plots below
  data(coastsCoarse)     # An alternative base map. Needs one fix:
  ind = which(coastsCoarse@lines[[94]]@Lines[[1]]@coords[,1] > 180)
  coastsCoarse@lines[[94]]@Lines[[1]]@coords[ind,1] = 180
  
  if(base.map=='countries')       base.map=countriesCoarse
  else  base.map = coastsCoarse
  
  # Handy constants
  n.bin = length(bins) - 1
  
  # Coordinate reference systems. Assuming WGS84 for GCD lat/lon. Some options for projection (besides default) 
  # are Mercator ("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") 
  # or Gall-Peters ("+proj=cea +lon_0=0 +lat_ts=45 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_def").
  # These "proj4" strings can be found online (e.g. http://www.spatialreference.org/)
  
  crs.unproj = CRS("+proj=longlat +datum=WGS84")
  crs.proj   = CRS(proj.str)
  
  
  # ----- Create composite
  if(class(PF)=="pfTransform") {
    cat("Creating composite...")
    # Run pfComposite. Not interested in the composite, but this will do the binning for us. 
    # (I've checked "by hand" and it is accurate and efficient.)
    COMP = pfComposite(PF, bins=bins, nboot=1, binning=T)
    
    if(is.null(id))
      id = pfSiteSel(ID=PF$params$IDn$SitesIDs)
    cat("done!\n\n")
  } 
  
  
  # ----- Data prep
  cat("Preparing data...")
  # Get lat/lon from GCD. The option to include them as inputs was added for quick testing.
  if(is.null(site.dat)) 
    site.dat = summary(id)
  
  # Convert data to spatial data into a SpatialPointsDataFrame (useful for doing the spatial stuff below)
  # Add binned + transformed charcoal data to original data frame. 
  dat = cbind(site.dat$LONGITUDE, site.dat$LATITUDE, t(COMP$BinnedData)) # A matrix
  dat = as.data.frame(dat) # A data frame
  names(dat) = c('lon','lat',paste('Z',1:n.bin, sep='')) # A data frame with named columns
  
  # Specify the columns that are coordinates, and then the coordinate ref. system. 
  # This converts the data frame to class 'SpatialPointsDataFrame', a special class in the 'sp' library. 
  coordinates(dat) = c('lon','lat')
  proj4string(dat) = crs.unproj
  cat("done!\n")
  
  # ----- Prediction grid. 
  cat("Defining grid...")
  # Get the extent of the full map
  full.extent.unproj = raster(nrows=1,ncols=1, crs=crs.unproj, ext=extent(grd.ext))
  
  if(proj.map) {   # Gall-peters (units are meters)
    # Get origin and dimension of projected map
    full.extent.proj = projectExtent(full.extent.unproj, crs.proj)
    
    grd.origin  = c(xmin(full.extent.proj), ymin(full.extent.proj)) # southeast corner
    grd.dim     = c(xmax(full.extent.proj)-xmin(full.extent.proj),
                    ymax(full.extent.proj)-ymin(full.extent.proj))
    
    if(is.null(grd.res))           # Use default of 1/100 of full extent
      grd.res     = grd.dim/100
    else if(length(grd.res)==1)    # Assume equal x/y resolution
      grd.res     = rep(grd.res,2)
    
    if(is.null(max.kr.dist)) 
      # Default maximum distance for kriging. Here, chosen to equal the distance from the center
      # to a corner of each grid cell.
      max.kr.dist = sqrt(sum(grd.res^2))/2
    
    # Project the original data and the coastlines (used later in plots)
    dat = spTransform(dat, crs.proj)
    base.map = spTransform(base.map, crs.proj)
    
  } else {
    # Unprojected (units are degrees, except max.kr.dist which is great-circle dist in km)
    
    # Get origin and dimension of unprojected map
    grd.origin  = c(xmin(full.extent.unproj), ymin(full.extent.unproj)) 
    grd.dim     = c(xmax(full.extent.unproj)-xmin(full.extent.unproj),
                    ymax(full.extent.unproj)-ymin(full.extent.unproj))
    
    if(is.null(grd.res))           # Use default of 5??
      grd.res     = c(5,5)
    else if(length(grd.res)==1)    # Assume equal x/y resolution
      grd.res     = rep(grd.res,2)
    
    if(is.null(max.kr.dist)) 
      # Default maximum distance for kriging. Here, a tricky choice since grid is defined in degrees 
      # but this parameter is in km. 
      max.kr.dist = 500
  }
  
  # Create prediction grid
  # Define "grid topology" (another 'sp' package thing)
  grd.top = GridTopology(grd.origin+grd.res/2, grd.res, ceiling(grd.dim/grd.res))
  
  if(proj.map) 
    grd = SpatialGrid( grd.top, crs.proj )
  else
    grd = SpatialGrid( grd.top, crs.unproj)
  cat("done!\n")
  
  
  # ----- Plotting
  cat("Plotting ", n.bin, " figures...", sep="")
  
  # Create plot directory if plotting
  if(!is.null(fig.base.name))
    dir.create(dirname(fig.base.name), recursive=T)
  
  # Loop over all bins...
  for(i in 1:n.bin) {
    cat(i,"...",sep="")
    
    # Open file connection, if saving the plot
    if(!is.null(fig.base.name))
      tiff(paste(fig.base.name, "_", bins[i], "-", bins[i+1], "bp.tif", sep=""), width=800, height=1000)
    
    # Calculate averages and variances within cells
    dat.mn = aggregate(dat[,i], grd, mean, na.rm=T)
    dat.var = aggregate(dat[,i], grd, function(x) sd(x)/length(x))
    
    
    # Convert to rasters (seems the simplest for testing, though probably not best)
    rast.mn  = as(dat.mn, "RasterLayer") # Convert to raster
    rast.var  = as(dat.var, "RasterLayer") # Convert to raster
    
    # Parameters
    par(mfrow=c(2,1), mar=c(4,4,2,3), bg='white') # subplot layout (mfrow = c(nrows,ncols)), margins, background color
    aspect = 1 # Aspect ratio (can make the map fit the figure size better)
    
    # Z (color scale) limits
    zlim.mn = c(-2,2) # can set to 'NULL' (no quotes) to default to range of data
    zlim.var = c(0,2)
    
    # Colors (uses an add-on library of color palettes). First argument is # of levels, then name of palette. 
    # rev() function switches its direction (the result of brewer.pal(...) is just a vector of color codes)
    
    cols.mn = rev(brewer.pal(11,'RdBu')) 
    cols.var = brewer.pal(9, 'YlOrRd')
    
    # The library has a convenience function to show all the palettes available:
    #           par(mfrow=c(1,1))
    #           display.brewer.all()
    
    # Make the plots
    plot(rast.mn,xaxs='i', legend=T, asp=aspect, col=cols.mn, zlim=zlim.mn, main=paste("Grid-cell mean: ", bins[i], "-", bins[i+1], " BP", sep="")) #, xlim=x.lim, ylim=y.lim )
    plot  (base.map, add=T, col=grey(0.5)) # Add world coastlines
    #         points(dat, pch=20, cex=0.25)
    
    plot(rast.var,xaxs='i', legend=T, asp=aspect, col=cols.var, zlim=zlim.var, main=paste("SE of the grid-cell mean: ", bins[i], "-", bins[i+1], " BP", sep=""))
    plot(base.map, add=T, col=grey(0.5))
    #         points(dat, pch=20, cex=0.25)
    
    #  Close file, if saving plot. 
    if(!is.null(fig.base.name)) dev.off()
  }
  cat("done!\n")
  
  
  
  
  cat("\nAll done!\n\n")
}


# ********************** Leftover for now...
# x.lim = c(-1.5e7, -1e7)             # X limits of map (in projection units)
# y.lim = c(5e6,1e7)                  # Y limits of map

# Remove NA rows (sites not spanning this bin); they cause trouble with the gstat functions and are of no use anyway. 
#         dat = dat[!is.na(dat$Z),]

# ----- Kriging
# Remove duplicate sampling locations. Required for kriging because otherwise we end up with a singular covariance matrix. (***a better solution would average values at replicate sampling locations, but this works for testing***)
#       zd = zerodist(dat)
#       dat = dat[-zd[,1],]

# Variogram fitting
# Calculate empirical variogram
#         dat.vgm = variogram(Z~1, dat)
#   
#       # Specify fit type and initial guesses
#         # This is going to fail for the projected data (units are way wrong), but that's OK since below we use the auto-fit version anyway. 
#           vgm.model = vgm(1, "Exp", 5000, 0.5)
#           vgm.fit = fit.variogram(dat.vgm, vgm.model)
#   
#       # Or fit using autofit
#         vgm.auto = autofitVariogram(Z~1, dat)$var_model
# 
#       # Plot and compare
#         par(mfrow=c(1,1))
#         x.lim.variogram = c(0,max(dat.vgm$dist)); y.lim=c(0,max(dat.vgm$gamma))
#         plot(dat.vgm$dist, dat.vgm$gamma, col=1, xlim=x.lim.variogram, ylim=y.lim)
#         lines(variogramLine(vgm.fit, max(dat.vgm$dist)), col=4)
#         lines(variogramLine(vgm.auto, max(dat.vgm$dist)), col=2)
#           legend("bottomright",col=c(1,4,2), lty=c(0,1,1), pch=c(1,NA,NA),c("empirical", "my fit", "auto fit"), bty='n')
# 
#     # Perform Kriging. Note high sensitivity to maxdist argument
#       dat.kr = krige( Z~1, locations=dat, newdata=grd, model=vgm.auto, maxdist=max.kr.dist)
#   
#     # Set cells with no sites to NA
#     if(!allow.interpolate) {
#       ind.nodata = which( !(1:nrow(dat.kr) %in% over(dat,grd)) )
#       dat.kr[ind.nodata,] = NA
#     }
# *************************************




