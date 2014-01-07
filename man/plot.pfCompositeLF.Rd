\name{plot.pfCompositeLF}
\alias{plot.pfCompositeLF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot.pfCompositeLF}
\description{
Plot pfCompositeLF object}
\usage{
\method{plot}{pfCompositeLF}(x, type = "ci", add = NULL, conf = c(0.05, 0.95), palette = "jet",
                 ...)
}
\arguments{
  \item{x}{
A "pfCompositeLF" object.}
  \item{type}{
Character, type of plot among "ci", "prctile", "density"}
\item{add}{
Character, add=NULL by default, add="sitenum" could be specified to plot the sites number in eah bin along with the composite curve.}
  \item{conf}{
Numeric, confidence levels.}
  \item{palette}{
Character, color palette used with type=c("prctile", "density") among "jet" and "BW".}
  \item{\dots}{\dots}
}

\author{
O. Blarquez}

\examples{
ID=pfSiteSel(Region=c("WNA0"),L12=1,Longlim=c(-160,-140))

TR=pfTransform(ID, method=c("MinMax","Box-Cox","MinMax","Z-Score"),
               BasePeriod=c(200,2000),QuantType="INFL")

COMP1=pfCompositeLF(TR, tarAge=seq(-50,4000,10), hw=200, nboot=100)

plot(COMP1)

## Note: comparing confidence intervals based on 100 replicates is not recommended
# (100 is used to decrease analysis time)

}
