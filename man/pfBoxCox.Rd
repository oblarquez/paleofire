\name{pfBoxCox}
\alias{pfBoxCox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Box-Cox transformation of Charcoal series}
\description{
Box-Cox transformation of charcoal series, the maximum likelihood estimation of lambda 
is derived from the boxcox.R function in the Venables and Ripley MASS library included in R 2.6.1 
}
\usage{
pfBoxCox(serie, alpha = 0.01, type = "BoxCox1964")
}
\arguments{
  \item{serie}{
A vector of charcoal values.}
  \item{alpha}{
Numeric, the "shift" parameter, default=0.01.}
  \item{type}{
Character, the Box-Cox transformation formulation, can be either "BoxCox1964" (default) for the original Box & Cox (1964) formulation, or "JohnDraper" for the John & Draper (1980) modulus transformation.}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{X}{Vector of transformed charcoal values}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Venables, W. N., Ripley, B. D., & Venables, W. N. (1994). Modern applied statistics with S-PLUS (Vol. 250). New York: Springer-verlag.

Box, G.E.P. & Cox, D. R.(1964) An analysis of transformations, Journal of the Royal Statistical Society, Series B, 26, 211-252.

John, J. A. & Draper N. R. (1980) Analternative family of transformations, Applied Statistics, 29, 190-197.
}
\author{
P. Bartlein}
\seealso{\code{\link{pfTransform}} }
\examples{
# Select a site
ID=pfSiteSel(SiteName="Pas-de-Fond")

# Extract data
A=pfExtract(ID)

B=pfBoxCox(A[,4],0.1)
plot(B,type="l")

}
