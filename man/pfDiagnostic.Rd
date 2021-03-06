% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pfDiagnostic.R
\name{pfDiagnostic}
\alias{pfDiagnostic}
\title{Print diagnostic pdf for individual transformed series}
\usage{
pfDiagnostic(
  ID,
  add = NULL,
  Age = 0,
  Interpolate = FALSE,
  method = "Box-Cox",
  BasePeriod = c(-100, 1e+09),
  span = 0.3,
  RunWidth = 500,
  RunQParam = 0.5,
  stlYears = 500,
  alpha = 0.01,
  type = "BoxCox1964",
  FileName = "Diagnostic.pdf",
  QuantType = "ALL"
)
}
\arguments{
\item{ID}{An object returned by \code{\link{pfSiteSel}} or
\code{\link{pfTransform}}}

\item{add}{An object returned by \code{\link{pfAddData}}}

\item{Age}{Numeric, if Interpolate=TRUE, Age is used to specified the ages
where the interpolation took place, If Age=0 the interpolated ages are
automatically specified using the median resolution of the record(s) If Age
is specified as a vector (e.g. Age=(from=0,to=10000, by=10)) the
interpolation took place at specified ages}

\item{Interpolate}{Logical, indicates wether data should be interpolated or
not, default=FALSE}

\item{method}{A character indicating the transformation method: "Z-Score",
Z-Score, "LOESS", Locally weighted regression, "SmoothSpline", Smoothing
spline, "Box-Cox", Box-Cox transformation, "MinMax", Minimax transformation,
"RunMed", Running median, "RunMean", Running mean, "RunQuantile", Running
quantile, "RunMin", Running min, "RunMax", Running max, "stl", Decompose a
time series into seasonal, trend and irregular components using loess, based
on \code{\link[stats]{stl}} function.}

\item{BasePeriod}{Numeric, a parameter specifying the base period for
calculating Z-score given in years BP (e.g. BasePeriod=c(0, 4000)), if empty
or unspecified the base period corresponds to record length.}

\item{span}{Numeric, the span parameter for the LOESS or Smoothing spline
methods}

\item{RunWidth}{Numeric, the width of the window for the"RunMed", "RunMean",
"RunQuantile", "RunMin", and "RunMax" methods in years.}

\item{RunQParam}{Numeric, the parameter specifying which quantile should be
calculated for the method "RunQuantile" (default=0.5 i.e. median).}

\item{stlYears}{Numeric, the bandwith for stl decomposition, default=500
years.}

\item{alpha}{Numeric, alpha value to add before BoxCox calculation, see
\code{\link{pfBoxCox}}.}

\item{type}{Character, the type of Box-Cox transformation, see
\code{\link{pfBoxCox}} for details}

\item{FileName}{Character, define output pdf file name e.g.
FileName="mydata.pdf"}

\item{QuantType}{Character, by default QuantType="INFL" and influx are
automatically calculated, otherwise use QuantType="NONE" (not recommended).}
}
\value{
\item{Filename.pdf }{A diagnostic file is printed, each sites being printed
on separate pages (specified using FileName="myfile.pdf"")}
}
\description{
Print diagnostic pdf for individual transformed series, successive
transformations could be specified (see example)
}
\examples{
\dontrun{
# Select boreal sites from Levavasseur 2012 PNV in Western North America
ID=pfSiteSel(continent=="North America", l12==1, long>=-160 & long<=-140)

# Print a diagnostic pdf for Box-Cox, Smoothed and Z-score tranformed data
# (base period = 200-2000 BP)
pfDiagnostic(ID,method=c("Box-Cox", "SmoothSpline","Z-Score"),
             span=0.3,BasePeriod=c(200,4000))

}
}
\author{
O. Blarquez
}
