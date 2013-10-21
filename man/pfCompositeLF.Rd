\name{pfCompositeLF}
\alias{pfCompositeLF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Produce a composite serie from multiple charcoal records using a local regression procedure (from the locfit package)}
\description{
Produce a composite serie from multiple charcoal local fits produced using the locfit procedure on bootstraped series subsamples, the sites charcoal values are prebinned prior to sites resampling. This procedure is equivalent to Daniau et al. 2012.
}
\usage{
pfCompositeLF(TR,hw=250,tarAge=NULL,binhw=NULL,nboot=1000,conf=c(0.05,0.95),
                       pseudodata=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TR}{
An object returned by pfTransform()}
  \item{tarAge}{
The target ages for prebinning given in years (e.g. tarAge=seq(0,10000,20)) 
If unspecified the sequence is defined as tarAge=seq(from=min age, to=max Age, by=median resolution) 
}
 \item{binhw}{
Bin half width for the prebinning procedure (use the same value as tarAge intervals for overlapping bins or tarAge intervals/2 for non-overlapping bins)}
  \item{nboot}{
A number specifying the number of bootstrap resamplings}
  \item{hw}{
The half window width for the locfit procedure (in years)}
\item{conf}{
Define confidence levels}
  \item{pseudodata}{
If TRUE 10 percent of the data is reflected at the top and the bottom of the resampled serie prior of each locfit regression in order to correct for the edge effect introduced by the local regression, see Cowling & Hall (1996).}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{out}{A matrix containing the c(0.05 0.5 0.95) quantiles issued from the bootstrap procedure and associated bin centres}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Cowling A, Hall P (1996) On pseudodata methods for removing boundary effects in kernel density estimation. Journal of the Royal Statistical Society, Series B 58(3):
551-563.}
\author{
O.Blarquez}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

ID=pfSiteSel(Latlim=c(30,70),Longlim=c(-100,-50))

TR=pfTransform(ID, method=c("MinMax","Box-Cox","MinMax","Z-Score"),
BasePeriod=c(200,2000),QuantType="INFL")

COMP1=pfCompositeLF(TR, tarAge=seq(-50,2000,10), hw=200, nboot=100)
COMP2=pfCompositeLF(TR, tarAge=seq(-50,2000,10), hw=200, nboot=100, pseudodata=TRUE)

dev.off()
par(mfrow=c(2,1))
plot(COMP1)
plot(COMP2)
## Note: comparing confidence intervals based on 100 replicates is not recommended
# (100 is used to decrease analysis time)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
