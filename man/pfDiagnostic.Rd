\name{pfDiagnostic}
\alias{pfDiagnostic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print diagnostic pdf for individual transformed series}
\description{
Print diagnostic pdf for individual transformed series, successive transformations could be specified (see example)}
\usage{
pfDiagnostic(IDn, add = NULL, Age = 0, Interpolate = FALSE, method
                 = "Box-Cox", BasePeriod = c(-100, 1e+09), span = 0.3,
                 RunWidth = 500, RunQParam = 0.5, stlYears = 500, alpha
                 = 0.01, type = "BoxCox1964", FileName =
                 "Diagnostic.pdf", QuantType = "ALL")
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
  \item{alpha}{
alpha value to add before BoxCox calculation}
  \item{type}{
The type of Box-Cox transformation, see pfBoxCox for details}
\item{FileName}{
Define output pdf file name e.g. FileName="mydata.pdf"}
\item{QuantType}{
To calculates influx set QuantType="INFL" otherwise set  QuantType="NONE"}


}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{Filename.pdf }{A diagnostic file is printed, each sites being printed on separate pages (specified using FileName="myfile.pdf"")}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
O. Blarquez}


\examples{
# Select boreal sites from Levavasseur 2012 PNV
ID=pfSiteSel(L12=1)

# Print a diagnostic pdf for Box-Cox, Smoothed and Z-score tranformed data 
# (base period = 200-2000 BP)
pfDiagnostic(ID,method=c("Box-Cox", "SmoothSpline","Z-Score"),
             span=0.3,BasePeriod=c(200,4000))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the

