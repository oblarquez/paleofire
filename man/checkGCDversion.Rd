\name{checkGCDversion}
\alias{checkGCDversion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check GCD package install}
\description{
Check if GCD package is installed and up to date during paleofire attach to ensure always using the most up to date GCD version. devtools package is required: on Windows install Rtools.exe depending on your R version http://cran.r-project.org/bin/windows/Rtools/
}
\usage{
checkGCDversion()
}
%- maybe also 'usage' for other objects documented here.
\details{
Last GCD database version is donwloaded and installed using:

\code{library(devtools)}

\code{install_github("GCD",username="paleofire",ref="master")}
}
\author{
O. Blarquez}

\examples{
checkGCDversion()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
