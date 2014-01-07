\name{plot.pfCircular}
\alias{plot.pfCircular}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot.pfCircular}
\description{
Plot circular block bootstrap percentiles.
}
\usage{
\method{plot}{pfCircular}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A "pfCircular" object.
}
\item{\dots}{\dots}

}
\author{
O. Blarquez}

\examples{
ID=pfSiteSel(Latlim=c(49,75),Longlim=c(6,50))
TR1=pfTransform(ID, method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,2000))

## Circular block bootstrapp
COMP=pfComposite(TR1, binning=TRUE, bins=seq(0,2000,100))
circ=pfCircular(COMP,conf=c(0.005,0.025,0.975,0.995),nboot=100)
plot(circ)
}
