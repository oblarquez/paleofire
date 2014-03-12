\name{pfTransform}
\alias{pfTransform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform charcoal data for unique to multiple series}
\description{
Charcoal data transformation, background estimation and homogenization for unique to multiple series, accepts objects returned by \code{\link{pfSiteSel}}.}
\usage{
pfTransform(ID = NULL, add = NULL, Interpolate = FALSE, Age = NULL, method
                 = "Z-Score", BasePeriod = c(-100, 1e+09), span = 0.3,
                 RunWidth = 500, RunQParam = 0.5, stlYears = 500, type
                 = "BoxCox1964", alpha = 0.01, QuantType = "INFL", 
                 MethodType = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ID}{
An object returned by \code{\link{pfSiteSel}} or \code{\link{pfTransform}}
}
\item{add}{
An object returned by \code{\link{pfAddData}}
}
  \item{Interpolate}{
Logical, idicates wether data should be interpolated or not, default=FALSE}
  \item{Age}{
Numeric, If Interpolate=TRUE, Age is used to specified the ages where the interpolation took place,
If Age=NULL (default) the interpolated ages are automatically specified using the median resolution of the record(s).
If Age is specified as a vector (e.g. Age=(from=0,to=10000, by=10)) the interpolation took place at specified ages.}
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
"stl", Decompose a time series into seasonal, trend and irregular components using loess, based on \code{\link[stats]{stl}} function.}
  \item{BasePeriod}{
Numeric, a parameter specifying the base period for calculating Z-score given in years BP (e.g. BasePeriod=c(0, 4000)), if empty or unspecified the base period corresponds to record length.}
  \item{span}{
Numeric, the span parameter for the LOESS or Smoothing spline methods}
  \item{RunWidth}{
Numeric, the width of the window for the"RunMed", "RunMean", "RunQuantile", "RunMin", and "RunMax" methods in years.}
  \item{RunQParam}{
Numeric, the parameter specifying which quantile should be calculated for the method "RunQuantile" (default=0.5 i.e. median).}
  \item{stlYears}{
Numeric, the bandwidth for stl decomposition, default=500 years.}
  \item{alpha}{
Numeric, alpha value to add before BoxCox calculation, see \code{\link{pfBoxCox}}.}
  \item{type}{
Character, the type of Box-Cox transformation, see \code{\link{pfBoxCox}} for details.}
\item{QuantType}{
Character, by default QuantType="INFL" and influx are automatically calculated, otherwise use QuantType="NONE" (not recommended).}
\item{MethodType}{
Character, by default (MethodType=NULL) imply that when for a specific site two charcoal unit exist the function pick the one define by pref_unit. By passing different arguments to MethodType user can modify the analysis to pick non preferred units by referring to more general methods for instance MethodType = "POLS" will choose charcoal records from pollen slides, or MethodType = "SIEV" sieved macro charcoal series. Type (paleofiredata); levels(paleofiredata$METHOD) for available methods.}
 \item{verbose}{
Logical, verbose or not...}
}

\value{
An object of the class "pfTransform".
}

\author{
O. Blarquez
}

\examples{
## Select the site Pas-de-Fond
ID=pfSiteSel(site_name=="Pas-de-Fond")

# Transform data sequentially using pfTransform function
tr=pfTransform(ID,method=c("MinMax","Box-Cox"))

## Plot transformed data for the first site
plot(tr$Age[,1],tr$TransData[,1],type="l")

}

