\name{plot.pfCompositeLF}
\alias{plot.pfCompositeLF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot.pfCompositeLF}
\description{
Plot pfCompositeLF object}
\usage{
\method{plot}{pfCompositeLF}(x, type = "ci", conf = c(0.05, 0.95), palette = "jet",
                 ...)
}
\arguments{
  \item{x}{
A "pfCompositeLF" object.}
  \item{type}{
Character, type of plot among "ci", "prctile", "density"}
  \item{conf}{
Numeric, confidence levels.}
  \item{palette}{
Character, color palette used with type=c("prctile", "density") among "jet" and "BW".}
  \item{\dots}{\dots}
}

\author{
O. Blarquez}

\examples{
ID=pfSiteSel(Latlim=c(30,70),Longlim=c(-100,-50))
TR1=pfTransform(ID, method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,2000))
COMP2=pfCompositeLF(TR1, tarAge=seq(-50,12000,10), hw=500,nboot=100)
plot(COMP2)

}
