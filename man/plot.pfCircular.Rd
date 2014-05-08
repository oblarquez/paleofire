\name{plot.pfCircular}
\alias{plot.pfCircular}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot.pfCircular}
\description{
Plot circular block bootstrap percentiles.
}
\usage{
\method{plot}{pfCircular}(x,ylim=NULL,xlim=NULL,ylab=NULL,xlab=NULL,main=NULL,text=FALSE,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A "pfCircular" object.
}
\item{ylim}{
Numeric, x axis limits.}
 \item{xlim}{
Numeric, y axis limits.}
\item{ylab}{
Character, y axis label.}
 \item{xlab}{
Character, x axis label.}
 \item{main}{
Character, title of the plot.}
 \item{text}{
Logical, text options.}
\item{\dots}{\dots}

}
\author{
O. Blarquez}

\examples{
ID=pfSiteSel(lat>49,lat<75,long>6,long<50)
TR1=pfTransform(ID, method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,2000))

## Circular block bootstrapp
COMP=pfComposite(TR1, binning=TRUE, bins=seq(0,2000,100))
circ=pfCircular(COMP,conf=c(0.005,0.025,0.975,0.995),nboot=100)
plot(circ)
}
