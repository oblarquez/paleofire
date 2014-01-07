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
An object of the class "pfSiteSel".}
  \item{type}{
Character, type of plot among "Map" or "Chronology".}
  \item{zoom}{
Character, zooming factor for type="Map": "Sites" or "World"}
  \item{\dots}{\dots}
}

\author{
O. Blarquez}

\examples{
ID=pfSiteSel(Region="ENA0")
plot(ID)

}

