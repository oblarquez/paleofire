\name{summary.pfSiteSel}
\alias{summary.pfSiteSel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
summary.pfSiteSel}
\description{
Return a summary table for an object of the class "pfSiteSel"}
\usage{
\method{summary}{pfSiteSel}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
 An object of the class "pfSiteSel".
 }
  \item{\dots}{
\dots}
}

\value{
Data.framce, returns the following informations: "ID_SITE", "LATITUDE", "LONGITUDE", "ELEV","MIN_EST_AGE", "MAX_EST_AGE", "NUM_DATING", "NUM_SAMP".
}

\author{
O. Blarquez}

\examples{
ID=pfSiteSel(ID=2)
summary(ID)
}

