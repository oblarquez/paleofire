\name{plot.pfComposite}
\alias{plot.pfComposite}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot.pfComposite}
\description{
Plot a pfComposite object}
\usage{
\method{plot}{pfComposite}(x, type = "ci", conf = c(0.05, 0.95), palette = "jet", smoothing = FALSE, spar = 0.3, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A pfComposite object}
  \item{type}{
Type of plot among "ci", "prctile", "density"}
  \item{conf}{
Confidence levels}
  \item{palette}{
Color palette used with type=c("prctile", "density")}
  \item{smoothing}{
If true smooth the resulting density from type="density"
}
  \item{spar}{
Smoothin parameter}
  \item{\dots}{
Other arguments passed to plot...
}
}

\author{
O. Blarquez}

\examples{
ID=pfSiteSel(Latlim=c(30,70),Longlim=c(-100,-50))
TR1=pfTransform(ID, method=c("MinMax","Box-Cox","ZScore"),BasePeriod=c(200,2000))
COMP=pfComposite(TR1, binning=TRUE, bins=seq(0,2000,20))
plot(COMP,type="density",smoothing=TRUE,spar=0.3)
}
