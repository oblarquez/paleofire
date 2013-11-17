\name{pfExtract}
\alias{pfExtract}
\title{
Extract charcoal data for a list of sites}
\description{
Extract charcoal data from an IDn object obtained by \code{\link{pfSiteSel}}}
\usage{
pfExtract(IDn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IDn}{
An object returned by \code{\link{pfSiteSel}}.}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{out}{A matrix of charcoal data with the following structure: out[,1]=Site identifiers,  out[,2]=Depths,  out[,3]=Estimated ages,  out[,4]=Charcoal data.}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
O. Blarquez}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## Retrieve a site
ID=pfSiteSel(SiteName="Pas-de-Fond")
## Or a group of sites (Western North America)
ID=pfSiteSel(Region=c("WNA0"))

## Extract data
A=pfExtract(ID)

# Plot the first site raw charcoal data
plot(A[A[,1]==ID$SitesIDS[1],3],A[A[,1]==ID$SitesIDS[1],4],type="l",main=ID$SiteNames[1],
     xlab="Age",ylab="raw Char")
}
