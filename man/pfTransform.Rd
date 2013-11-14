\name{pfTransform}
\alias{pfTransform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform charcoal data for unique to multiple series}
\description{
Transform charcoal data for unique to multiple series, accepts objects returned by pfSiteSel() and pfTransform for successive data transformation}
\usage{
pfTransform(IDn, add = NULL, Interpolate = FALSE, Age = 0, method
                 = "Z-Score", BasePeriod = c(-100, 1e+09), span = 0.3,
                 RunWidth = 500, RunQParam = 0.5, stlYears = 500, type
                 = "BoxCox1964", alpha = 0.01, QuantType = "ALL")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IDn}{
An object returned by pfSiteSel() or pfTransform()
}
\item{add}{
Add data returned by pfAddData
}
  \item{Interpolate}{
Logical, idicates wether data should be interpolated or not}
  \item{Age}{
If Interpolate=TRUE, Age is used to specified the ages where the interpolation took place,
If Age=0 the interpolated ages are automatically specified using the median resolution of the record(s)
If Age is specified as a vector (e.g. Age=(0,10000, by=10)) the interpolation took place at specified ages}
  \item{method}{
A character indicating the transformation method:
"Z-Score", Z-Score, 
"LOESS", Locally weighted regression,
"SmoothSpline", Smoothing spline,
"Box-Cox", Box-Cox transformation,
"MinMax", Minimax transformation,
"RunMed", Running median,
"RunMean", Running mean,
"RunQuantile", Running quantile,
"RunMin", Running min,
"RunMax", Running max,
"stl", Decompose a time series into seasonal, trend and irregular components using loess, based on stl function.}
  \item{BasePeriod}{
A parameter specifying the base period for calculating Z-score given in years BP (e.g. BasePeriod=c(0, 4000)), if empty or unspecified the base period corresponds to record length}
  \item{span}{
The span parameter for the LOESS or Smoothing spline methods}
  \item{RunWidth}{
The width of the window for the"RunMed", "RunMean", "RunQuantile", "RunMin", and "RunMax" methods in years}
  \item{RunQParam}{
The parameter specifying which quantile should be calculated for the method "RunQuantile" 
(default=0.5 i.e. median)}
  \item{stlYears}{
The bandwith for stl decomposition, default=500 years}
  \item{type}{
The type of Box-Cox transformation, see pfBoxCox for details}
\item{alpha}{
alpha value to add before BoxCox calculation}
\item{QuantType}{
To calculates influx set QuantType="INFL" otherwise set  QuantType="NONE"}

}

\value{
An object of the class "pfTransform""
}

\author{
O. Blarquez
}

\examples{
## Select sites from the temperate deciduous/broadleaf forest
# (PNV 3 and 5 in Ramankutty and Foley 1999)
ID=pfSiteSel(RF99=c(3,5))

# Transform data sequentially using pfTransform function
tr=pfTransform(ID,method=c("MinMax","Box-Cox"))

## Plot transformed data for the first site
plot(tr$Age[,1],tr$TransData[,1],type="l")

}

