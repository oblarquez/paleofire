\name{pfSiteSel}
\alias{pfSiteSel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GCD sites selection methods}
\description{
Main function used for site selection, uses data strored in data(paleofiresites) to perform site selection according to multiple criterion, those criterions could be either geographic, based on series attributes (e.g. # of datings), or on sites attributes (e.g. biome)}
\usage{
pfSiteSel(ID = NULL, Latlim = c(-360, 360), Longlim = c(-360,
                 360), Biome = NULL, DateInt = NULL, Country = NULL,
                 Region = NULL, SiteName = NULL, PrefUnit = NULL,
                 Elevation = c(-1e+05, 1e+05), QuantType = NULL, L12=NULL, RF99=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{ID}{
Specify sites identifiers if known (ID_SITE in paleofiresites)}
  \item{Latlim}{
A vector specifying the latitudes limits
}
  \item{Longlim}{
A vector specifying the longitudes limits
}
  \item{Biome}{
A vector specifying the biome in wich the sites should be part
}
  \item{DateInt}{
The mean dating interval (see above example)}
  \item{Country}{
Country code}
  \item{Region}{
Region code}
  \item{SiteName}{
The site name}
  \item{PrefUnit}{
Charcoal units}
  \item{Elevation}{
Site elevation}
  \item{QuantType}{
Define the quantity type measure e.g. QuantType="INFL" for influx}
  \item{RF99}{Potential Natural Vegetation following Ramankutty and Foley (1999), see ?paleofiresites for details}
  \item{L12}{Potential Natural Vegetation following Levavasseur et al. (2012), see ?paleofiresites for details}
    }
    
\details{
The criterion values could be inspected by accessing data(paleofiresites)}
\value{
An object of the class pfSiteSel
}

\author{
O. Blarquez}

\examples{
library(paleofire)
require(GCD)

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
