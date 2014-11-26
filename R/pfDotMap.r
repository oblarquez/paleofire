#  v08
#  v06

pfDotMap = function(TR, bins, 
                    fig.base.name=NULL, base.map='coasts',
                    grd.res=5, grd.ext=c(-180,180,-90,90), 
                    proj4="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", n.boot=1000,
                    cx.minsize=0.3, cx.mult=1
) {
  if (!requireNamespace("rworldmap", quietly = TRUE)) {
    install.packages("rworldmap")
  }
  require("rworldmap")
  
  
  # ---------------- TEST BLOCK
  # Easier to test without running the code as a function. Comment everything above here (function definition) and 
  # uncomment this code to run the whole thing as a standard script. The load() line needs to point 
  # to where a pre-made pfTransform object has been saved. You will also need to load the two additional 
  # functions at the end of this file (just run them once at the beginning of each session). 
  # And technically, don't run the '}' that closes the main function definition (though I think if you do it will 
  # run everything and just give a harmless error at the end.
#   # 
#   rm(list=ls())
#   
#   #load('~/Drive/GPWG_MapPaper/Krige Maps/All_GPCD_Transformed_v2.rdata')
#   TR =res3
#   bins              = seq(0,2000,1000)
#   fig.base.name     = '/Users/Olivier/Desktop/'
#   base.map          = 'coasts'
#   grd.res           = 5
#   grd.ext           = c(-180,180,-90,90)
#   proj4             = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "  # Corrsponds to unprojected maps
#   n.boot            = 10    # too small, but OK for a teest
#   cx.minsize        = 0.3   # minimum dot size
#   cx.mult           = 1     # multiplicative factor for scaling all dots
#   
  # ---------------- END TEST BLOCK
  
  
  # ----- Libraries
  
  
  # ----- Load base map
  countriesCoarse<-coastsCoarse<-NULL
  rm(countriesCoarse);rm(coastsCoarse)
  
  data(countriesCoarse,envir = environment())  # A dataset in rworldmap used in the plots below
  data(coastsCoarse,envir = environment())     # An alternative base map. Needs one fix:
  ind = which(coastsCoarse@lines[[94]]@Lines[[1]]@coords[,1] > 180)
  coastsCoarse@lines[[94]]@Lines[[1]]@coords[ind,1] = 180
  
  # Select and transform base map based on inputs (only two options currently)
  if(base.map=='countries') {
    base.map = countriesCoarse
  } else {
    base.map = coastsCoarse
  }
  base.map = sp::spTransform(base.map, sp::CRS(proj4))
  
  
  # ----- Create composite
  if(class(TR)=="pfTransform") {
    cat("Creating composite...")
    # Run pfComposite. Not interested in the composite, but this will do the binning for us. 
    # (I've checked "by hand" and it is accurate and efficient.)
    COMP = pfComposite(TR, bins=bins, nboot=1000, binning=T)
    CHAR = t(COMP$BinnedData)
    n.bin = length(bins) - 1
    cat("done!\n")
  } else {
    stop("Input PF must be a pfTransform object.")
  }
  
  
  # ----- Get lat/lon from GCD
  cat("Retrieveing site coordinates from GCD...")
  # A little awkward, but seems good to rely on existing paleofire functions for retrieving site coordinates,
  # and as far as I can tell this is the way to do it at present.

  # Lookup info from the database
  paleofiresites=NULL; rm(paleofiresites)  
  data(paleofiresites,envir = environment())
  site.dat=paleofiresites[paleofiresites$id_site %in% TR$params$ID$id_site, ]
  
  
  # Extract coordinates, including in radians for use in the distance function
  sites.lon = site.dat$long
  sites.lonrad = sites.lon*pi/180
  sites.lat = site.dat$lat
  sites.latrad = sites.lat*pi/180
  
  # Define n.site
  n.site = nrow(site.dat)
  
  cat("done!\n")
  
  
  # ----- Define prediction grid
  if(length(grd.res)==1)    # Assume equal x/y resolution if single number given
    grd.res     = rep(grd.res,2)
  
  # Find grid cell centers
  grd.lon = seq( grd.ext[1]+grd.res[1]/2, grd.ext[2], grd.res[1])
  grd.lat = seq( grd.ext[3]+grd.res[2]/2, grd.ext[4], grd.res[2])
  
  # Expand to obtain every combination of lon/lat
  grd.lonlat = expand.grid(grd.lon,grd.lat)
  names(grd.lonlat) = c("lon", "lat") # column names for convenience
  n.grd = nrow(grd.lonlat)
  
  # Again, will need lat/lon in radians later  
  grd.lonlat.rad = grd.lonlat*pi/180
  
  
  # ------ Calculate distances from sites to grid cell centers
  cat("\nCalculating distances...\n")
  # Output space and progress bar definition
  dists = matrix(NA, nrow=n.grd, ncol=n.site)
  pb = txtProgressBar(0,n.grd,style=2)
  
  # Loop over every grid cell, save distance from the grid cell center to every site
  for(i in 1:n.grd) {   # i=1
    dists[i,] = haverdist(grd.lonlat.rad$lon[i], grd.lonlat.rad$lat[i],sites.lonrad,sites.latrad)
    setTxtProgressBar(pb, i)
  } 
  close(pb) # Close progress bar
  
  
  # ----- Compute stats for each grid cell
  cat("\nComputing stats for each grid cell X time slice...\n")
  #  As discussed at AGU (Marlon, Bartlein, Higuera, Kelly), a simple way is to search within a radius 
  # around each grid cell center, equal to the greatest distance from a grid cell center to its most distant corner. 
  # This should occur at the equator, where grid cells are largest. Note that this conservative approach will
  # result in many sites falling within multiple grid boxes--even at the equator, 
  # the defined circles will overlap near the edges of the grid boxes. At higher latitudes,
  # the grid cells are much smaller, so overlap will be considerably greater. 
  # There are alternatives, like using a grid that is irregular in terms of lat/lon, 
  # or changing the area of grid cells depending on latitude. But all have their tradeoffs (we thought), 
  # and this one is simple. 
  
  # Find the max distance just discussed. The center-to-corner distance for a cell at the origin should work:
  max.dist = haverdist(0,0, (grd.res[1]/2)*(pi/180), (grd.res[2]/2)*(pi/180))
  
  # Set aside space for some grid-cell statistics
  grd.n   = matrix(0, nrow=n.grd, ncol=n.bin)
  grd.mean  = matrix(NA, nrow=n.grd, ncol=n.bin)
  grd.lCI = matrix(NA, nrow=n.grd, ncol=n.bin)
  grd.uCI = matrix(NA, nrow=n.grd, ncol=n.bin)
  dat.n.contributions = matrix(0, nrow=n.site, ncol=n.bin) 
  
  # Now loop over each time bin X grid cell and compute desired stats.
  pb = txtProgressBar(0,n.grd*n.bin,style=2) # progress bar
  for(j in 1:n.bin) {  # i=1887; j=1
    for(i in 1:n.grd) {
      ind = which(dists[i,]<max.dist)         # all sites within range
      ind = ind[which(!is.na(CHAR[ind,j]))] # remove NA (sites that don't span bin j)
      if(length(ind)>0) {
        grd.n[i,j]    = length(ind)            # number of sites contributing
        grd.mean[i,j] = mean(CHAR[ind,j])      # mean CHAR
        CI = bootCI(CHAR[ind,j], nboot=n.boot) # CI
        grd.lCI[i,j] = CI[1] # split into lower/upper to keep variables 2-D
        grd.uCI[i,j] = CI[2]
        
        # For each site that contributed to this grid cell, increment dat.n.contributions
        dat.n.contributions[ind,j] = dat.n.contributions[ind,j] + 1
      }
      setTxtProgressBar(pb, (j-1)*n.grd+i) # update progress bar
    }
  }
  close(pb) # close progress bar

  grd.site.ind = list()
  for(i in 1:n.grd) {
      grd.site.ind[[i]] = which(dists[i,]<max.dist)         # all sites within range
  }


# ---------- Plot
  cat("\nPlotting ", n.bin, " figures...\n", sep="")
  pb = txtProgressBar(0,n.bin,style=2) # progress bar
  plotlist = spgrdlist = spsitelist = list()

  
  # ----- Loop over all bins, one plot for each
  for(j in 1:n.bin) {  #     j=1
    # ----- Setup
    #       cat(j,"...",sep="")
    setTxtProgressBar(pb, j) # update progress bar
    
    # Open file connection, if saving the plot
    if(!is.null(fig.base.name))
      pdf(paste(fig.base.name, "_", bins[j], "-", bins[j+1], "bp.pdf", sep=""), 
          width=12, height=12)
    
    # Convert stats to spatial data frames. Might not be necessary but sure makes it easy to use spplot() below.
    sp.grd = cbind(grd.lonlat, grd.n[,j], grd.mean[,j], grd.lCI[,j], grd.uCI[,j])
    names(sp.grd) = c("lon","lat","sitesPerCell","mean.CHAR","CI.lower","CI.upper")
    sp::coordinates(sp.grd) = c('lon','lat')
    sp::proj4string(sp.grd) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "
    
    sp.sites = data.frame(sites.lon, sites.lat, dat.n.contributions[,j])
    names(sp.sites) = c("lon","lat","cellsPerSite")
    sp::coordinates(sp.sites) = c('lon','lat')
    sp::proj4string(sp.sites) = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "
    
    # Now project the data to be plotted
    sp.grd     = sp::spTransform(sp.grd, sp::CRS(proj4))
    sp.sites   = sp::spTransform(sp.sites, sp::CRS(proj4))
    
    # Get x/y lims (could add option to override later). Only purpose currently is that sites vs. grid have 
    # different bounding boxes, which is the default extent for spplot(). I.e. without setting all the same, 
    # the bottom-right plot will have different extent than other two.
    x.lim = sp::bbox(sp.grd)[1,]
    y.lim = sp::bbox(sp.grd)[2,]

    # ----- Create mean plot
      # Define colors and cut locations 
#       cols = c("#CA0020","#F4A582","#F7F7F7","#92C5DE","#0571B0") # from colorbrewer
#       cols = c("#CA0020","#F4A582",grey(0.9),"#92C5DE","#0571B0") # modified from colorbrewer
      cols = c("#0571B0","#92C5DE",grey(0.9),"#F4A582","#CA0020") # modified from colorbrewer
      cuts = seq(-2.5,2.5,by=1) # Defines range and resolution of color scale
        
      # Determine cuts for sizing point.
        # Specify symbol sizes for the two classes
        cx.sizes = cx.mult*c(0.5,1) 
        
        # Assign symbol size based on whether CI contain 0 
        cx = ifelse(sp.grd$CI.lower>0 | sp.grd$CI.upper<0, max(cx.sizes), min(cx.sizes))
        
        # The previous line will produce NA for cells with n=1 since CI are undefined. Give these "non-significant" symbol size by default.
        cx[which(sp.grd$sitesPerCell==1)] = min(cx.sizes)
          
      # Create plot object (actually plotted later)
      mean.plot = 
        sp::spplot(sp.grd, 'mean.CHAR', xlim=x.lim, ylim=y.lim,
          cuts=cuts, colorkey=T, col.regions=cols, cex=cx, edge.col=grey(0.3), lwd=0.3,
          scales=list(draw=T), sp.layout=list("sp.lines",base.map,col=grey(0.8)),
          main=paste("Charcoal Influx z-Scores: ", bins[j], "-", bins[j+1], " BP", sep="")) 


    # ----- Plot Number of sites per grid cell
      # Generate dot sizes/colors and corresponding key. 
        # Specify scale and legend. Would be good to automate this but it's pretty tricky to produce a good general algorithm. So, for now, hard-coding symbol sizes / labels that work well for global map at 5?? resolution. 
        cuts      = c(0,1,5,10,20,1000) # Where to divide symbol sizes
        cols      = grey(0.2)   # Can be replaced by a vector if different colors are desired
        cx.legend = c("1", "2-5", "6-10", "11-20",">20") # legend text
          n.cx = length(cuts)-1   # number of bins represented
      
        # Define sizes of data points and legend entries by dividing data into bins and scale to range [cx.minsize,1]*cx.mult. This is a pretty good range for symbol 'cex' sizes, although can be modified with the cx.mult and cx.minsize arguments. 'cut(..., labels=F)' returns integer classes from 1:n.cx. These are the same sizes to use in the key. 
        cx.key = ( ((1:n.cx)-1)*(cx.mult-cx.minsize)/(n.cx-1) + cx.minsize )
        cx = cx.key[ cut(sp.grd$sitesPerCell, cuts, labels=F) ]


      # Adjust scale so that the low end corresponds to specified minimum symbol size
      ind.non0 = which(cx>0) # Don't want to change size 0 (== not plotted)
      cx.key = cx.key + cx.minsize - min(cx[ind.non0]) 
      cx[ind.non0] = cx[ind.non0] + cx.minsize - min(cx[ind.non0])


      # Create plot object (actually plotted later)
      sitesPerCell.plot = 
        sp::spplot(sp.grd, 'sitesPerCell', xlim=x.lim, ylim=y.lim, scales=list(draw=F), 
          cex=cx, cex.key=cx.key, legendEntries=cx.legend, cuts=cuts, 
          col.regions=cols, edge.col="white", lwd=0.3,
          sp.layout=list("sp.lines",base.map,col=grey(0.8)), key.space="right",
          main="Number of sites per grid cell")

 
    # ----- Plot Number of grid cells contributed per site
      # Generate dot sizes/colors and corresponding key. 
        # Specify scale and legend. Would be good to automate this but it's pretty tricky to produce a good general algorithm. So, for now, hard-coding symbol sizes / labels that work well for global map at 5?? resolution. 
        cuts   = c(0,1,2,3,4,100) # Where to divide symbol sizes
        cols      = grey(0.2)   # Can be replaced by a vector if different colors are desired
        cx.legend = c("1", "2", "3", "4",">4") # legend text
          n.cx = length(cuts)-1   # number of bins represented
      
        # Define sizes of data points and legend entries. 
        cx.key = ( ((1:n.cx)-1)*(cx.mult-cx.minsize)/(n.cx-1) + cx.minsize )
        cx = cx.key[ cut(sp.sites$cellsPerSite, cuts, labels=F) ]

      # Adjust scale so that the low end corresponds to specified minimum symbol size
      ind.non0 = which(cx>0) # Don't want to change size 0 (== not plotted)
      cx.key = cx.key + cx.minsize - min(cx[ind.non0]) 
      cx[ind.non0] = cx[ind.non0] + cx.minsize - min(cx[ind.non0])
        
      # Create plot object (actually plotted later)
      cellsPerSite.plot = 
        sp::spplot(sp.sites, 'cellsPerSite', xlim=x.lim, ylim=y.lim,
          cex=cx, cex.key=cx.key, legendEntries=cx.legend, cuts=cuts, 
          scales=list(draw=F), col.regions=cols, edge.col="white", lwd=0.3,
          sp.layout=list("sp.lines",base.map,col=grey(0.8)), key.space="right",
          main="Number of grid cells influenced by each site")

  
    
    # ----- Create time series plot
      timeSeries.dat = data.frame(
                age  = rep(bins,each=2)[2:(2*n.bin+1)],
                char = rep(COMP$Result$MEAN, each=2),
                lCI  = rep(COMP$Result[,3], each=2),
                uCI  = rep(COMP$Result[,4], each=2) )


      timeSeries.plot =       
        xyplot( char~age, data=timeSeries.dat, 
          ylim=c(-0.05,0.05)*diff(range(timeSeries.dat[,3:4]))+range(timeSeries.dat[,3:4]),
          panel = function(x,y, ...) {
            ind.j = (2*j-1):(2*j)
            x.j   = x[ind.j]
            y.j   = y[ind.j]
            lCI.j = timeSeries.dat$lCI[ind.j]
            uCI.j = timeSeries.dat$uCI[ind.j]
            
            panel.polygon(c(x,rev(x)), c(timeSeries.dat$lCI, rev(timeSeries.dat$uCI)), 
              col=grey(0.8), border=FALSE)
            panel.polygon(c(x.j,rev(x.j)), c(lCI.j, rev(uCI.j)), 
              col=rgb(1,0,0,0.5), border=FALSE)

            panel.xyplot(x,y, type='l', col=1)
            panel.xyplot(x.j,y.j, type='l', col=2)
          })


    # ----- Produce plots
      print(mean.plot, position=c(0,0.5,1,1), more=T)
      print(timeSeries.plot, position=c(0,0.3,1,0.5), more=T)
      print(sitesPerCell.plot, position=c(0,0,0.5,0.3), more=T)
      print(cellsPerSite.plot, position=c(0.5,0,1,0.3))


    # Close connection to external figure
    if(!is.null(fig.base.name))  dev.off()

    # ----- Store plot objects
      plotlist[[j]] = list("mean"=mean.plot,"sitesPerCell"=sitesPerCell.plot,"cellsPerSite"=cellsPerSite.plot,"timeSeries"=timeSeries.plot)
      spgrdlist[[j]] = sp.grd
      spsitelist[[j]] = sp.sites
      

  } # End loop over all bins

  # ----- Return
    output = list(COMP=COMP, bins=bins, sp.grd=spgrdlist, sp.sites=spsitelist, grd.site.ind=grd.site.ind, site.dat=site.dat, plots=plotlist)
    return(output)
  
cat("\nAll done!\n\n")
  

} # End main function definition



# ---------- Auxilliary functions
# Haversine distance function (http://en.wikipedia.org/wiki/Haversine_formula)
# R is the approximate radius of the earth, defaulting to 6371 km
# Note: at least on pair of inputs (i.e., either (lon1,lat1) or (lon2/lat2)) must represent a single location. 
# The other may be vectors representing many locations, but the function does not work to find all pairwise 
# distances between two sets (n>1) of locations. It wouldn't be difficult to code the pairwise version, 
# but I think it wouldn would require a loop so it wouldn't save any run time anyway.

haverdist = function(lon1, lat1, lon2, lat2, R=6371) {
  return( 
    2*R*asin(sqrt( (sin((lat2-lat1)/2))^2 + cos(lat1)*cos(lat2)*((sin((lon2-lon1)/2))^2) ))
  )
}

# Bootstrap CI function.
# Actually should work for various functions and CI probs but this code currently only relies 
# on the defaults (mean, 95% CI).
bootCI = function(x, nboot=1000, fun=mean, ..., probs=c(0.025,0.975)) {
  nx = length(x)
  
  if(nx<2) {
    return( rep(NA, length(probs)) )
  } else {
    xboot = matrix(x[sample(1:nx, nboot*nx, replace=T)], nrow=nboot, ncol=nx)
    return( quantile(apply(xboot, 1, fun, ...), probs=probs, na.rm=T) )
  }
}
