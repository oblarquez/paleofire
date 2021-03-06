% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contiguous.R
\name{plot.contiguous}
\alias{plot.contiguous}
\title{Plot "contiguous" object}
\usage{
\method{plot}{contiguous}(x, ylim = NULL, xlim = NULL, ...)
}
\arguments{
\item{x}{An object returned by contiguous}

\item{ylim}{Numeric, ylim for the graph}

\item{xlim}{Numeric, xlim for the graph}

\item{\dots}{\dots{}}
}
\value{
A plot.
}
\description{
Plot an object returned by contiguous, plot contiguous cores (or sites) in
green (T) and non-contiguous cores in red (F).
}
\examples{
\dontrun{
x=pfSiteSel(lat>12,lat<60,long<(-50),long>-140)
cont=contiguous(x)
plot(cont)
}
}
\seealso{
\code{\link{contiguous}}
}
\author{
O. Blarquez
}
