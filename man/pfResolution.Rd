\name{pfResolution}
\alias{pfResolution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates age resolution indicators for charcoal records}
\description{
Calculates age resolution indicators for charcoal records selected using \code{\link{pfSiteSel}} or \code{\link{pfInteractive}} functions.}
\usage{
pfResolution(ID, AgeLim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ID}{
An object of the class "pfSiteSel"}
  \item{AgeLim}{
Numeric, defines age limits for age resolution calculations (e.g. AgeLim=c(-50,6000))}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{data.frame}{A data frame with the following informations: ID_SITE, SITE_NAME, Median Resolution of the record, Mean Resolution and Standard seviation}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
O. Blarquez}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
ID1=pfSiteSel(Latlim=c(45,90), Longlim=c(-100,-50))
Res=pfResolution(ID1,AgeLim=c(-50,8000))
head(Res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
