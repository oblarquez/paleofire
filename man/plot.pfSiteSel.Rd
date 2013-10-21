\name{plot.pfSiteSel}
\alias{plot.pfSiteSel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot.pfSiteSel
}
\description{
Plot an object of the class "pfSiteSel"}
\usage{
\method{plot}{pfSiteSel}(x, type = "Map", zoom = "Sites", \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of the class "pfSiteSel"}
  \item{type}{
"Map" or "Chronology"}
  \item{zoom}{
Zooming factor for type="Map": "Sites" or "World"}
  \item{\dots}{
Other arguments to be passed to plot}
}

\author{
O. Blarquez}

\examples{
library(paleofire)
require(GCD)

ID=pfSiteSel(ID=c(1,3,13,67))
plot(ID)
plot(ID,type="Chronology")


}

