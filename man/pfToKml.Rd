\name{pfToKml}
\alias{pfToKml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Export selected site to Google Earth kml format}
\description{
Export sites selected using pfSiteSel function to Google Earth kml format.
}
\usage{
pfToKml(x, file = "NULL")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of the class "pfSiteSel"}
  \item{file}{
File location and name with kml extension e.g. file="/Users/Olivier/Desktop/truc.kml"}
}

\value{
No value returned.
}

\author{
O. Blarquez}

\examples{
\dontrun{
x=pfSiteSel(id_site==222)
pfToKml(x,file="site222.kml")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
