\name{plot.pfKruskal}
\alias{plot.pfKruskal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a "pfKruskal" object.}
\description{
Plot a "pfKruskal" object using boxplots and showing significant differences between the periods using letters.}
\usage{
\method{plot}{pfKruskal}(x, trend = FALSE, outliers = FALSE,  xlim = NULL, ylim = NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object returned by \code{\link{pfKruskal}}.}
  \item{trend}{
Logical, show trend using linear regression?}
  \item{outliers}{
Logical, show outliers?}
  \item{xlim}{
Numeric, x axis limits.}
  \item{ylim}{
Numeric, y axis limits.}
  \item{\dots}{
...}
}
\value{
Return a ggplot2 "gg" object.
}
\details{
If two periods share the same letter their rank (median) is not significantly different at the confidence level specified by \code{alpha}. If not, equality could be rejected at the confidence level specified by \code{alpha}.
}
\author{
O. Blarquez}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pfKruskal}}
}
\examples{
## Composite charcoal record for Western Boreal North America:
ID=pfSiteSel(id_region=="WNA0" & l12==1)
plot(ID)
## Transform data
res3=pfTransform(ID,method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,4000))

## Composite
comp=pfComposite(res3,bins=seq(from=-500,to=12500,by=1000))
plot(comp)

## Kruskal Wallis Anova
comparison=pfKruskal(comp)

plot(comparison)

# p=plot(comparison)
# require(ggplot2)
# p+ggtitle("my title")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
