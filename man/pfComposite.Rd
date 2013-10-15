\name{pfComposite}
\alias{pfComposite}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Produce a composite serie from multiple charcoal records}
\description{
Produce a composite serie from multiple charcoal records using bootstrap resampling,
the sites charcoal values are binned and the mean in each bin is calculated prior the bootstrap procedure. This procedure is equivalent to Power et al. 2008.
}
\usage{
pfComposite(TR, bins = NULL, nboot = 1000, binning = TRUE, conf =
                 c(0.05, 0.95))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TR}{
An object returned by pfTransform()}
  \item{bins}{
The sequence for binning given in years (e.g. bins=seq(0,10000,200)) 
If unspecified the sequence is defined as bins=seq(from=min age, to=max Age, by=median resolution) 
}
  \item{nboot}{
A number specifying the number of bootstrap resamplings}
\item{binning}{
Set to TRUE (default) for binning, if transformed data are first interpolated this agument can be set to FALSE (no binning)}
\item{conf}{
Define confidence levels}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{out}{A matrix containing the c(0.05 0.5 0.95) quantiles issued from the bootstrap procedure and associated bin centres}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
O.Blarquez}

\examples{

## Composite charcoal record for North America:
ID=pfSiteSel(Region=c("ENA0","WNA0"))

## Transform data
res3=pfTransform(ID,method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,4000),Interpolate=FALSE)

## Composite
comp=pfComposite(res3,bins=seq(0,20000,200))
plot(comp)
}
