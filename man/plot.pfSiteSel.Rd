\name{plot.pfSiteSel}
\alias{plot.pfSiteSel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot.pfSiteSel
}
\description{
Plot an object of the class "pfSiteSel"}
\usage{
\method{plot}{pfSiteSel}(x, add = NULL, type = "Map", zoom = "Sites", pch =
                 "|", xlim = NULL, ylim = NULL, cex = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of the class "pfSiteSel".}
\item{add}{
An object returned by pfAddData (optional).}
  \item{type}{
Character, type of plot among "Map" or "Chronology".}
  \item{zoom}{
Character, zooming factor for type="Map": "Sites" or "World"}
\item{pch}{
Pointer type see \code{\link[graphics]{plot}}.}
\item{xlim}{
Numeric, x axis limits.}
 \item{ylim}{
Numeric, y axis limits.}
\item{cex}{
Numeric, size of points.}
  \item{\dots}{\dots}
}

\author{
O. Blarquez}

\examples{
ID=pfSiteSel(id_region=="ENA0", long>-100)
plot(ID,zoom="world")

}

