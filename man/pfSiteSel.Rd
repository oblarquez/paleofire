\name{pfSiteSel}
\alias{pfSiteSel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GCD sites selection methods}
\description{
Main function used for site selection, uses data strored in data(paleofiresites) to perform site selection according to multiple criterion, those criterions could be either geographic, based on series attributes (e.g. # of datings), or on sites attributes (e.g. biome).}
\usage{
pfSiteSel(ID = NULL, Latlim = c(-360, 360), Longlim = c(-360,
                 360), Biome = NULL, DateInt = NULL, Country = NULL,
                 Region = NULL, SiteName = NULL, PrefUnit = NULL,
                 Elevation = c(-1e+05, 1e+05), QuantType = NULL, L12=NULL, RF99=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{ID}{
Numeric, specify sites identifiers if known (ID_SITE in paleofiresites).}
  \item{Latlim}{
Numeric, a vector specifying the latitudes limits.
}
  \item{Longlim}{
Numeric, a vector specifying the longitudes limits.
}
  \item{Biome}{
Numeric, a vector specifying the biome in which the sites should be located.
}
  \item{DateInt}{
Numeric, the mean dating interval (see example below).}
  \item{Country}{
Factor, country code.}
  \item{Region}{
Factor, region code.}
  \item{SiteName}{
Character, the sites names.}
  \item{PrefUnit}{
Factor, charcoal units.}
  \item{Elevation}{
Numeric, a vector specifying sites elevations limits.}
  \item{QuantType}{
Character, define the quantity type measure e.g. QuantType="INFL" for influx.}
  \item{RF99}{Numeric, Potential Natural Vegetation following Ramankutty and Foley (1999), see ?paleofiresites for details.}
  \item{L12}{Numeric, Potential Natural Vegetation following Levavasseur et al. (2012), see ?paleofiresites for details.}
    }
    
\details{
The criterion values could be inspected by accessing data(paleofiresites) in GCD package.}
\value{
An object of the class pfSiteSel
}
\seealso{\code{\link[GCD]{paleofiresites}}}
\author{
O. Blarquez}

\examples{
## Sites selection examples

## Select all sites
ID=pfSiteSel()

## Site in the Biome #8
ID=pfSiteSel(Biome=8)

## Sites in North America by geographic location
ID=pfSiteSel(Latlim=c(25,75),Longlim=c(-150,-45))
## By region criterion
ID=pfSiteSel(Region=c("ENA0","WNA0"))

## Pas-de-Fond site
ID=pfSiteSel(SiteName="Pas-de-Fond")

## Sites with on average one dating point every 250 yrs
ID=pfSiteSel(DateInt=250)

## Sites between 0, 100 m elevation in Asia
ID=pfSiteSel(Elevation=c(0,100),Region="ASIA")

## Explore individual criterions
data(paleofiresites)
unique(paleofiresites$ID_REGION)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
