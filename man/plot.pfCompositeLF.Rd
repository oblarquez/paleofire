\name{plot.pfCompositeLF}
\alias{plot.pfCompositeLF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot.pfCompositeLF}
\description{
Plot pfCompositeLF object}
\usage{
\method{plot}{pfCompositeLF}(x, type = "ci", \dots)
}
\arguments{
  \item{x}{
pfCompositeLF object}
  \item{type}{
Simple "ci" is implemented (more to come)}
  \item{\dots}{
Other arguments passed to plot...}
}

\author{
O. Blarquez}

\examples{
ID=pfSiteSel(Latlim=c(30,70),Longlim=c(-100,-50))
TR1=pfTransform(ID, method=c("MinMax","Box-Cox","ZScore"),BasePeriod=c(200,2000))
COMP2=pfCompositeLF(TR1, bins=seq(-50,12000,20), hw=500,nboot=100)
plot(COMP2)

}
