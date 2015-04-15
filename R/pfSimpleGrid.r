#' Produce simple gridded maps of paleofire data
#' 
#' Produce gridded map graphics representing spatial variability in charcoal 
#' data from the Global Charcoal Database.
#' 
#' Takes any pfTransform object as input, and allows any set of one or more
#' time bins to be specified for plotting (one plot per bin).
#' 
#' Results will be plotted on a regular lon/lat grid. To determine which sites
#' contribute to each grid cell value, the code searches within a specified
#' great circle distance (i.e. on the surface of the globe) around each grid
#' cell center. To avoid missing any sites, the distance is set equal to the
#' greatest distance from a grid cell center to its most distant corner, which
#' occurs at the equator where grid cells are largest. This conservative
#' approach will result in many sites falling within multiple grid boxes. At
#' all latitudes, the defined radii will overlap near the edges of the grid
#' boxes. At higher latitudes, the lon/lat grid cells are physically much
#' smaller, so overlap will be considerably greater. There are alternatives,
#' like using a grid that is irregular in terms of lon/laton, or changing the
#' area of grid cells depending on latitude. But all have their tradeoffs, and
#' this one is simple.
#' 
#' Current version produces plots of mean CHAR, number of sites per grid cell,
#' and number of grid cells contributed to by each site (due to overlapping
#' radii described above). The mean plot additionally shows points in two
#' sizes, representing those mean values whose 95"\%" confidence intervals do
#' (small dots) or do not (large dots) contain zero. Finally, a time series is
#' plotted in each figure with the current time bin highlighted.
#' 
#' @param TR An object returned by \code{\link{pfTransform}}
#' @param tarAge Numeric, the target ages for prebinning given in years (e.g.
#' tarAge = seq(0, 10000, 20)). If unspecified the sequence is defined as
#' tarAge=seq(from=min age, to=max Age, by=median resolution).
#' @param hw Numeric, the half window width for the locfit procedure (in
#' years).
#' @param binhw Numeric, bin half width for the prebinning procedure (use the
#' same value as tarAge intervals for overlapping bins or tarAge intervals/2
#' for non-overlapping bins).
#' @param fig.base.name Character sequence representing the base name for the
#' figures. Can be preceded by a path as long as all directories in the path
#' exist. One figure will be produced for each time bin, with years (and file
#' suffix) appended to the base name automatically. A value of \code{NULL}
#' (default) causes figures to be plotted to the current device in sequence.
#' @param grd.res,grd.ext Desired grid resolution and extent in degrees. If
#' \code{grd.res} is a single number, the grid will be defined with equal
#' lon/lat resolution; a two-element vector (lon,lat) can also be supplied for
#' unequal resolution. \code{grd.ext} is specified as a vector of the form
#' \code{c(min-lon,max-lon,min-lat,max-lat)}.
#' @param grd.lonlat A data frame of coordinates for every grid cell center, to
#' be used in cases where an irregular grid is desired. Columns must be named
#' 'lon' and 'lat'. If specified, grd.res and grd.ext are ignored. Note that
#' this option could have undesirable results for unusual grid definitions. In
#' particular, the maximum radius for including sites in a grid cell is always
#' calculated at the equator. For a regular lon/lat grid, this guarantees all
#' sites will be included in at least one cell, because equatorial cells are
#' largest at the equator. If an irregular grid is specified such that this is
#' not true, the maximum radius calculated could lead to sites excluded from
#' all cells. In this case a warning is printed but the function proceeds
#' anyway.
#' @param base.map Currently, either \code{'coasts'} or \code{'countries'} to
#' choose which base map (from required library \code{'rworldmap'}) to be
#' plotted as the base map for all plots. Could easily be modified to accept
#' any SpatialPolygons object.
#' @param proj4 proj.4 string representing the desired projection for plotted
#' maps. Default is unprojected. See \url{http://www.spatialreference.org} to
#' look up the string for your favorite projections.
#' @param n.boot Number of bootstrap replicates to use when creating confidence
#' intervals around each grid-cell mean. In each time bin X grid cell
#' combination, replicates consist of composite z-score values for that bin,
#' randomly sampled (with replacement) from sites within the grid cell (see
#' 'Details' for precise description of sites included in each cell). I.e., no
#' temporal bootstrapping is done here, so that bootstrap CI reflect only
#' spatial variability.
#' @param cx.minsize,cx.mult Parameters that crudely adjust plotted dot size.
#' cx.minsize defines the minimum cex applied to any point in any map, cx.mult
#' scales all points by an equivalent factor.
#' @return Plots are produced on the current device or in pdf files defined by
#' \code{fig.base.name}. In addition, a named list of useful objects is
#' returned:
#' 
#' \item{COMP}{ The binned composite generated for plotting.  } \item{bins}{
#' The list of bin endpoints.  } \item{sp.grd}{ A
#' \code{\link[sp]{SpatialPointsDataFrame-class}} object containing all the
#' grid-level statistics produced and plotted (mean influx value, bootstrap
#' confidence interval, and number of sites per grid cell).  } \item{sp.sites}{
#' A \code{\link[sp]{SpatialPointsDataFrame-class}} object representing the
#' number of grid cells influenced by each site.  } \item{plots}{ A list with
#' one element for each bin. These elements are themselves named lists of
#' trellis objects representing each of the plots produced ("mean",
#' "sitesPerCell", "cellsPerSite", "timeSeries"). Note that these objects can
#' be edited to some degree with the \code{\link[lattice]{update.trellis}}
#' function, and plotted or used in layouts as any other trellis graphics can.
#' }
#' @author R. Kelly
#' @references Power, M., J. Marlon, N. Ortiz, P. Bartlein, S. Harrison, F.
#' Mayle, A. Ballouche, R. Bradshaw, C. Carcaillet, C. Cordova, S. Mooney, P.
#' Moreno, I. Prentice, K. Thonicke, W. Tinner, C. Whitlock, Y. Zhang, Y. Zhao,
#' A. Ali, R. Anderson, R. Beer, H. Behling, C. Briles, K. Brown, A. Brunelle,
#' M. Bush, P. Camill, G. Chu, J. Clark, D. Colombaroli, S. Connor, A. L.
#' Daniau, M. Daniels, J. Dodson, E. Doughty, M. Edwards, W. Finsinger, D.
#' Foster, J. Frechette, M. J. Gaillard, D. Gavin, E. Gobet, S. Haberle, D.
#' Hallett, P. Higuera, G. Hope, S. Horn, J. Inoue, P. Kaltenrieder, L.
#' Kennedy, Z. Kong, C. Larsen, C. Long, J. Lynch, E. Lynch, M. McGlone, S.
#' Meeks, S. Mensing, G. Meyer, T. Minckley, J. Mohr, D. Nelson, J. New, R.
#' Newnham, R. Noti, W. Oswald, J. Pierce, P. Richard, C. Rowe, M. Sanchez
#' Goni, B. Shuman, H. Takahara, J. Toney, C. Turney, D. Urrego-Sanchez, C.
#' Umbanhowar, M. Vandergoes, B. Vanniere, E. Vescovi, M. Walsh, X. Wang, N.
#' Williams, J. Wilmshurst, and J. Zhang. 2008. Changes in fire regimes since
#' the Last Glacial Maximum: an assessment based on a global synthesis and
#' analysis of charcoal data. Climate Dynamics 30:887-907.
#' @examples
#' 
#' \dontrun{
#' ## Composite charcoal record for North America:
#' ID=pfSiteSel(id_region==c("WNA0"), l12==1 & long<(-130))
#' plot(ID)
#' 
#' ## Transform data
#' res3=pfTransform(ID,method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,4000))
#' 
#' ## Plot maps for 1000-yr bins spanning 3-0 kBP
#' # dev.new(width=10,height=10) # A big plot area helps. 
#' gridmap = pfSimpleGrid( TR=res3, tarAge=seq(0,2000,1000), hw=500, ext=c(-170,-80,40,80))
#' summary(dotmap)
#' 
#' # Plot the mean map from the first time bin
#' # newmap = update(dotmap$plots[[1]]$mean, main="A relabeled map")
#' # newmap
#' }
#' 

pfSimpleGrid = function(TR, tarAge, hw, binhw=0.5*mean(diff(tarAge)), fun=mean,
                    n.boot=0, prob.CI=c(0.025,0.975), test.val=0, 
                    proj4="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",
                    res=5, ext=c(-180,180,-90,90), 
                    fig.file.name=NULL, show.plots=TRUE, title.text="",
                    cols = NULL, cuts = NULL, zlim=NULL,
                    base.map='coasts', base.map.col=grey(0.7), base.map.lwd=0.5
                    
) {
 
# ---------------- TEST BLOCK
# 
# # Easier to test without running the code as a function. Comment everything above here (function definition) and 
# # uncomment this code to run the whole thing as a standard script. The load() line needs to point 
# # to where a pre-made pfTransform object has been saved. You will also need to load the two additional 
# # functions at the end of this file (just run them once at the beginning of each session). 
# # And technically, don't run the '}' that closes the main function definition (though I think if you do it will 
# # run everything and just give a harmless error at the end.
# rm(list=ls())
# 
# load('/Work/Research/GPWG/pfDotMap working/All_GPCD_Transformed_v2.rdata')
# bins              = seq(-500,21500,1000)
# show.plots        = T
# fig.file.name     = '/Work/Research/GPWG/pfDotMap working/SimpleGrid Maps vTest7.pdf' # NULL for no saved plot
# base.map          = 'coasts'
# res           = 5
# ext           = c(-180,180,-90,90)
# proj4         = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "
# n.boot            = 100    # too small, but OK for a teest
# zlim              = c(-2.5,2.5)
# fun               = mean
# cols              = c("#0571B0","#92C5DE",grey(0.9),"#F4A582","#CA0020")  # Mean
# cuts = seq(-2.5,2.5,by=1) # Defines range and resolution of color scale
# 
# # zlim              = NULL
# # fun               = "count"
# # cols              = brewer.pal(9, 'YlOrRd')  # others
# # cuts              = NULL
# 
# prob.CI           = c(0.025,0.975)
# test.val          = 0
# 
# base.map.col      = grey(0.7)
# base.map.lwd      = 0.5
# title.text        = ""
# ---------------- END TEST BLOCK
  
  # ----- Constants
  # proj4 string for the GCD data
    proj4.GCD = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "

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
    # No need for bootstraps--we're only going to use the actual composite for now. Will bootstrap by site later. 
    COMP = pfCompositeLF(TR, tarAge=tarAge, hw=hw, binhw=binhw, nboot=1)
    CHAR = t(COMP$BinnedData)
    n.bin = length(tarAge)
    cat("done!\n")
  } else {
    stop("Input PF must be a pfTransform object.")
  }
  
  
  # ----- Get site locations from GCD
  cat("Retrieveing site coordinates from GCD...")
  # Lookup info from the database
    paleofiresites=NULL; rm(paleofiresites)  
    data(paleofiresites,envir = environment())
    site.dat=paleofiresites[paleofiresites$id_site %in% TR$params$ID$id_site, ]

  # Convert data to spatial data into a SpatialPointsDataFrame
    # Add binned + transformed charcoal data to original data frame. 
      dat = data.frame(lon=site.dat$long, lat=site.dat$lat, CHAR)

  # Specify the columns that are coordinates, and then the coordinate ref. system. This converts the data frame to class 'SpatialPointsDataFrame', a class in the 'sp' library. 
    sp::coordinates(dat) = c('lon','lat')
    sp::proj4string(dat) = proj4.GCD
    
  # NA locations
    na.mat = is.na(CHAR)
  cat("done!\n")
  
  
  # ----- Define prediction grid
  grd = raster(extent(ext))     # Extent
    res(grd) = res              # Resolution
    proj4string(grd) = proj4    # Projection


  # ----- Calculate desired stat for each time bin
  sg.rast = rasterize(dat, grd, field=names(dat), fun=fun, na.rm=T)
    # If don't specify "names(dat)" then get an extra layer "ID". 

    # Convert to stack (will be layer if only one time bin)
    sg.rast = stack(sg.rast)

  # ----- Bootstrap CI
  if(n.boot>0) {
    # Bootstrap CI function.
    # Actually should work for various functions and CI probs but this code currently only relies 
    # on the defaults (mean, 95% CI).
    bootTest = function(x, ...) {
      x = na.omit(x)
      nx = length(x)
  
      if(nx<2) {
        return( NA )
      } else {
        xboot = matrix(x[sample(1:nx, n.boot*nx, replace=T)], nrow=n.boot, ncol=nx)
        quants = quantile(apply(xboot, 1, fun, ...), probs=prob.CI)
        return(ifelse( quants[1]>max(test.val) | quants[2]<min(test.val), 1, 0))
#         return(quants)
      }
    }
  
    isSingle = function(x, ...) ifelse(length(na.omit(x))==1,1,NA)
#     boot.lCI = function(x,...) bootTest(x, prob=min(prob.CI))
#     boot.uCI = function(x,...) bootTest(x, prob=max(prob.CI))

    # Even though it doesn't need to be done in a loop, do so in order to print progress
    sg.rast.lCI = sg.rast.uCI = sg.rast.1 = sg.rast.sig = NA*sg.rast
    timer = Sys.time()
    for(i in 1:n.bin) {
      cat(i,"...")
#       sg.rast.lCI[[i]] = rasterize(dat, grd, field=names(dat)[i], fun=boot.lCI)
#       sg.rast.uCI[[i]] = rasterize(dat, grd, field=names(dat)[i], fun=boot.uCI)
#       sg.rast.sig[[i]] = (sg.rast.lCI[[i]]>=max(test.val)) | (sg.rast.uCI[[i]]<=min(test.val))
      sg.rast.sig[[i]] = rasterize(dat, grd, field=names(dat)[i], fun=bootTest)
      sg.rast.1[[i]] = rasterize(dat, grd, field=names(dat)[i], fun=isSingle)
      print(summary(as.factor(values(sg.rast.sig[[i]]))))
    }
    sg.rast.sig[ sg.rast.sig<1 ] = NA
    Sys.time()-timer
  }


  # ----- Handle default colors
  if(is.null(zlim)) zlim = range(values(sg.rast), na.rm=T)
  if(length(cuts)==1) cuts = seq(min(zlim), max(zlim), length=cuts+1)
  
  if(is.null(cols) & is.null(cuts)) {
    cuts = seq(min(zlim), max(zlim), length=10)
    cols = rev(topo.colors(length(cuts)-1))
  } else if(is.null(cols)) {
    cols = rev(topo.colors(length(cuts)-1))
  } else if(is.null(cuts)) {
    cuts = seq(min(zlim), max(zlim), length=length(cols)+1)
  }


  # ----- Make Plots
  cat("\nPlotting ", n.bin, " figures...", sep="")
  
  # Can call spplots on the whole raster brick, which creates a separate panel for each layer (in this case, one per time bin). E.g.

    # sg.plots = spplot(sg.rast, at=cuts, col.regions=cols, scales=list(draw=T),layout=c(1,1))
    # 
    # # Then can update the "proper" way for trellis objects. E.g.:
    #   sg.plots = update(sg.plots, 
    #               strip=strip.custom(factor.levels=
    #                 paste0(title.text, bins[1:n.bin], "-", bins[2:(n.bin+1)], " BP")),
    #               sp.layout=list("sp.lines",base.map,col=base.map.col,lwd=base.map.lwd))
      
  # Pretty snazzy. But a little confusing if like me you're not very familiar with lattice graphics. Also, not the way we set up pfDotMaps. 
  
  # So for now, making a separate trellis plot for each time bin, storing them in a list. 
  sg.plots = list()
  base.map.layout = list("sp.lines",base.map,col=base.map.col,lwd=base.map.lwd)
  pb = txtProgressBar(0,n.bin,style=2) # progress bar
  for(j in 1:n.bin) {  #     j=1
    setTxtProgressBar(pb, j) # update progress bar

    site.dots.layout = list("sp.points", dat[ !na.mat[,j],j], col=1, pch=16, cex=0.3)

    sg.plots[[j]] = sp::spplot(sg.rast[[j]], at=cuts, col.regions=cols, scales=list(draw=T),
      main = as.list(paste0(title.text, bins[j], "-", bins[j+1], " BP")))
    
    # Should be better way, but can't figure out how to add multiple extra layers in separate steps. This works...
    if(n.boot>0) {
      if(all(is.na(values(sg.rast.sig[[j]])))) {
        sg.sig.layout = list("sp.text", "(no sig vals)", loc=c(0,0))
      } else {
        sg.sig.layout = list("sp.polygons", rasterToPolygons(sg.rast.sig[[j]]), first=F)
      }
      
      if(all(is.na(values(sg.rast.1[[j]])))) {
        sg.1.layout = list("sp.text", "(no single-site cells)", loc=c(0,0))
      } else {
        sg.1.layout = list("sp.polygons", rasterToPolygons(sg.rast.1[[j]]), col=2, first=F)
      }
      sg.plots[[j]] = update(sg.plots[[j]], 
        sp.layout=list(base.map.layout,sg.1.layout,sg.sig.layout,site.dots.layout))
    } else {
      sg.plots[[j]] = update(sg.plots[[j]], sp.layout=list(base.map.layout, site.dots.layout))
    }
  }


  # ----- Print plots to screen and/or file
  # Plot to screen
    if(show.plots) print(sg.plots)

  # Save to file
    if(!is.null(fig.file.name)) {
      pdf(fig.file.name, width=12, height=12)
      print(sg.plots)
      dev.off()
    }
  cat("done!\n")


  # ----- Return
    output = list(COMP=COMP, bins=bins, sg.rast=sg.rast, sg.plots=sg.plots, site.dat=site.dat)
    return(output)
  
cat("\nAll done!\n\n")
} # End main function definition

