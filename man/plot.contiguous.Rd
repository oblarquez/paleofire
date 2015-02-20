\name{plot.contiguous}
\alias{plot.contiguous}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot "contiguous" object}
\description{
Plot an object returned by contiguous, plot contiguous cores (or sites) in green (T) and non-contiguous cores in red (F). }
\usage{
\method{plot}{contiguous}(x, ylim = NULL, xlim = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object returned by contiguous}
  \item{ylim}{
Numeric, ylim for the graph}
  \item{xlim}{
Numeric, xlim for the graph}
  \item{\dots}{
\dots
}
}

\value{
A plot.
}

\author{
O. Blarquez}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{contiguous}}
}
\examples{
x=pfSiteSel(lat>12,lat<60,long<(-50),long>-140)
cont=contiguous(x)
plot(cont)
}
