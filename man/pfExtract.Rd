\name{pfExtract}
\alias{pfExtract}
\title{
Extract charcoal data for a list of sites}
\description{
Extract charcoal data from an object returned by \code{\link{pfSiteSel}}}
\usage{
pfExtract(ID)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ID}{
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
ID=pfSiteSel(site_name=="Pas-de-Fond")
## Or a group of sites (Western North America)
ID=pfSiteSel(id_region==c("WNA0"))

## Extract data
A=pfExtract(ID)

# Plot the first site raw charcoal data
plot(A[A[,1]==ID$id_site[1],3],A[A[,1]==ID$id_site[1],4],type="l",main=ID$site_name[1],
     xlab="Age",ylab="raw Char")
}
