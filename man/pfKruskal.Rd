\name{pfKruskal}
\alias{pfKruskal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Analyse composite records by a Kruskal-Wallis ANOVA}
\description{
The function applies a Kruskal-Wallis ANOVA on binned data issued from a "pfComposite" object (of directly on "pfTransform" objects), in order to test the difference in biomass burning activity between different time periods.}
\usage{
pfKruskal(data, p.adj = "none", alpha = 0.05, bins=NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
An object returned by \code{\link{pfComposite}} or \code{\link{pfTransform}}.}
  \item{p.adj}{
Method for adjusting p values (see \code{\link[stats]{p.adjust}}). Includes: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr" and "none" (default).}
  \item{alpha}{
Numeric, confidence level.}
\item{bins}{
Numeric, bins to use if a "pfTransform" object is provided.}
  \item{verbose}{
Logical, verbose or not...}
}

\value{
A "pfKruskal" object containing multiple comparison results.
}

\author{
O. Blarquez}

\seealso{
\code{\link{plot.pfKruskal}},\code{\link[agricolae]{kruskal}}
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
}

