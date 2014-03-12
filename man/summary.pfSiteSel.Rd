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
Data.frame, returns the following informations: "id_site", "lat", "long"
"elev", "min_est_age", "max_est_age", "num_dating",  "date_int", "num_samp", 
"l12", "rf99".
}

\author{
O. Blarquez}

\examples{
ID=pfSiteSel(id_site==2)
summary(ID)
}

