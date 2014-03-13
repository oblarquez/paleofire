\name{pfMinMax}
\alias{pfMinMax}
\title{
MiniMax transformation of a charcoal serie}
\description{
MiniMax transformation of a charcoal serie}
\usage{
pfMinMax(serie)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{serie}{
Numeric, a vector of charcoal values.}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{out}{A vector of minimax transformed values.}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\seealso{\code{\link{pfTransform}}}
\author{
O. Blarquez}

\examples{
## Retrieve a site
ID=pfSiteSel(site_name=="Pas-de-Fond")
## Or a group of sites (Western North America)
ID=pfSiteSel(id_region==c("WNA0"))

## Extract data
A=pfExtract(ID)

## Plot the first site raw charcoal data
par(mfrow=c(1,2))
plot(A[A[,1]==ID$id_site[1],3],A[A[,1]==ID$id_site[1],4],type="l",main=ID$site_name[1],
     xlab="Age",ylab="raw Char")
## Minimax transformation
B=pfMinMax(A[A[,1]==ID$id_site[1],4])
## Plot the first site Minimax transformed charcoal data
plot(A[A[,1]==ID$id_site[1],3],B,type="l",main=ID$site_name[1],
     xlab="Age",ylab="Minimax")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

