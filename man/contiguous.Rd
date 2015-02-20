\name{contiguous}
\alias{contiguous}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Are cores sampled contiguously?}
\description{
The function checks wether cores have been sampled contiguously or with a depth resolution <1cm.}
\usage{
contiguous(x,threshold=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of the class "pfSiteSel"}
  \item{threshold}{
Numeric, threshold for considering two samples as contiguous (default=1cm)}
}

\value{
Summary table of sites with the added contiguous logical column (TRUE--FALSE)
}

\author{
O. Blarquez}

\seealso{
\code{\link{pfResolution}}
}
\examples{
x=pfSiteSel(lat>12,lat<60,long<(-50),long>-140)
contiguous(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
