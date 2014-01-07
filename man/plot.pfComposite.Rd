\name{plot.pfComposite}
\alias{plot.pfComposite}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot.pfComposite}
\description{
Plot a pfComposite object}
\usage{
\method{plot}{pfComposite}(x, type = "ci", conf = c(0.05, 0.95), palette = "jet",
                 add = "NONE", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A "pfComposite" object.}
  \item{type}{
Character, type of plot among "ci", "prctile", "density"}
  \item{conf}{
Numeric, confidence levels.}
  \item{palette}{
Character, color palette used with type=c("prctile", "density") among "jet" and "BW".}
 \item{add}{
Character, add="NONE" by default, add="sitenum" could be specified to plot the sites number in eah bin along with the composite curve.}
\item{\dots}{\dots}
}

\author{
O. Blarquez}

\examples{
## Composite charcoal record for North America:
ID=pfSiteSel(Region=c("WNA0"),L12=1)
## Transform data
res3=pfTransform(ID,method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,4000),Interpolate=FALSE)

## Composite
comp=pfComposite(res3,bins=seq(0,20000,200))
plot(comp,type="density",smoothing=TRUE,spar=0.3)
}
