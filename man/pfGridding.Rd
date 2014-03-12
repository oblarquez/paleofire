\name{pfGridding}
\alias{pfGridding}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Produce gridded maps of transformed charcoal values.}
\description{
The function uses weighted spatio-temporal interpolation to produce gridded maps of transformed charcoal values. Spatial grids are used to interpolate transformed charcoal values for a key period defined by \code{Age}. For each grid cell the function search charcoal sites located in a radius defined by \code{distance_buffer} from the grid centre and at an elevation within a range defined by \code{elevation_buffer} from the mean elevation of the cell. Then the function search  for charcoal samples within a temporal range from the key date defined by \code{time_buffer}. Finally a tricube distance weighting function is applied to each sample by considering it spatial distance to the grid centre and it temporal distance to the key date. This approach that weight samples according to their spatio-temporal location also down-weight charcoal sites that are poorly sampled.
}
\usage{
pfGridding(data, cell_sizex = NULL, cell_sizey = NULL, age = 0, 
           cell_size = NULL, time_buffer = NULL, 
           distance_buffer = NULL, raster_extent = NULL, 
           elevation_buffer = NULL, proj4 = NULL, 
           sea_mask = FALSE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
An object returned by \code{\link{pfTransform}}.}
  \item{cell_sizex}{
Numeric, grid cell width (m).}
  \item{cell_sizey}{
Numeric, grid cell height (m).}
  \item{age}{
Numeric, key date (years BP).}
  \item{cell_size}{
Numeric, grid cell size (bypass cell_sizex and cell_sizey and produce square cells).}
  \item{time_buffer}{
Numeric, temporal distance (years) from the key date to search for charcoal samples.}
  \item{distance_buffer}{
Numeric, spatial distance from the grid centres to search for charcoal samples (m).}
  \item{raster_extent}{
Numeric, define custom extent for the analysis such as raster_extent = c(xmin, xmax, ymin, ymax)}
  \item{elevation_buffer}{
Numeric, elevation range from the mean grid cell elevation to search for charcoal sites.}
  \item{proj4}{
String, proj.4 string representing the desired projection for plotted maps. Default is Robinson ("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"). See \url{http://www.spatialreference.org} to look up the string for your favorite projections.
}
  \item{sea_mask}{
Logical, mask cells falling in the sea.}
  \item{verbose}{
Logical, verbose or not...}
}

\value{
A "pfGridding" object (list) that could be plotted using \code{\link{plot.pfGridding}}. 
}

\author{
O.Blarquez}

\seealso{
\code{\link{plot.pfGridding}}, \code{\link{pfTransform}}, \code{\link{pfDotMap}}}
\examples{
ID=pfSiteSel(id_region=="ENA0", l12==1, long>-85)

TR=pfTransform(ID,method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,4000))

p=pfGridding(TR,age=1000)
summary(p)

# require(raster)
# plot(p$raster)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
