\name{paleofire-package}
\alias{paleofire-package}
\alias{paleofire}
\docType{package}
\title{
paleofire: A package for the Global Charcoal Database 
}
\description{
The paleofire package provides tools to extract and analyse charcoal sedimentary data stored in the Global Charcoal Database. Main functionalities includes data extraction and sites selection, transformation and interpolation of the charcoal records as well as compositing.
}
\details{
\tabular{ll}{
Package: \tab paleofire\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2013-05-23\cr
License: \tab GPL (>=2)\cr
}
}
\section{Author(s)}{Global Paleofire Working Group <paleofire@gmail.com>}
\section{Maintainer}{Olivier Blarquez <blarquez@gmail.com>}

\references{
Daniau, A. L., P. J. Bartlein, S. P. Harrison, I. C. Prentice, S. Brewer, P. Friedlingstein, T. I. Harrison-Prentice, J. Inoue, K. Izumi, J. R. Marlon, S. Mooney, M. J. Power, J. Stevenson, W. Tinner, Andri, M., J. Atanassova, H. Behling, M. Black, O. Blarquez, K. J. Brown, C. Carcaillet, E. A. Colhoun, D. Colombaroli, B. A. S. Davis, D. D'Costa, J. Dodson, L. Dupont, Z. Eshetu, D. G. Gavin, A. Genries, S. Haberle, D. J. Hallett, G. Hope, S. P. Horn, T. G. Kassa, F. Katamura, L. M. Kennedy, P. Kershaw, S. Krivonogov, C. Long, D. Magri, E. Marinova, G. M. McKenzie, P. I. Moreno, P. Moss, F. H. Neumann, E. Norstrom, C. Paitre, D. Rius, N. Roberts, G. S. Robinson, N. Sasaki, L. Scott, H. Takahara, V. Terwilliger, F. Thevenon, R. Turner, V. G. Valsecchi, B. Vanniere, M. Walsh, N. Williams, and Y. Zhang. 2012. Predictability of biomass burning in response to climate changes. Global Biogeochem. Cycles 26:GB4007.

Power, M., J. Marlon, N. Ortiz, P. Bartlein, S. Harrison, F. Mayle, A. Ballouche, R. Bradshaw, C. Carcaillet, C. Cordova, S. Mooney, P. Moreno, I. Prentice, K. Thonicke, W. Tinner, C. Whitlock, Y. Zhang, Y. Zhao, A. Ali, R. Anderson, R. Beer, H. Behling, C. Briles, K. Brown, A. Brunelle, M. Bush, P. Camill, G. Chu, J. Clark, D. Colombaroli, S. Connor, A. L. Daniau, M. Daniels, J. Dodson, E. Doughty, M. Edwards, W. Finsinger, D. Foster, J. Frechette, M. J. Gaillard, D. Gavin, E. Gobet, S. Haberle, D. Hallett, P. Higuera, G. Hope, S. Horn, J. Inoue, P. Kaltenrieder, L. Kennedy, Z. Kong, C. Larsen, C. Long, J. Lynch, E. Lynch, M. McGlone, S. Meeks, S. Mensing, G. Meyer, T. Minckley, J. Mohr, D. Nelson, J. New, R. Newnham, R. Noti, W. Oswald, J. Pierce, P. Richard, C. Rowe, M. Sanchez Goni, B. Shuman, H. Takahara, J. Toney, C. Turney, D. Urrego-Sanchez, C. Umbanhowar, M. Vandergoes, B. Vanniere, E. Vescovi, M. Walsh, X. Wang, N. Williams, J. Wilmshurst, and J. Zhang. 2008. Changes in fire regimes since the Last Glacial Maximum: an assessment based on a global synthesis and analysis of charcoal data. Climate Dynamics 30:887-907.

}

\keyword{ fire, global, charcoal, sediments, paleo}
\seealso{
\url{http://gpwg.org}
}
\examples{
rm(list=ls())
rm(list=ls())

library(paleofire)

## Interactive sites selection:
# ID=pfInteractive()

## Site selection using criterions
# DateInt parameter is used to set the mean interval which is required between two
# dating points (ex 14C) for sites to be selected for a complete list of criterions
# that can be used see pfSiteSel function

ID=pfSiteSel(Latlim=c(30,70),Longlim=c(-100,-50), DateInt=3000)
plot(ID,zoom="sites")

## Filter sites based on sample number using summary function
sumID=summary(ID)
sites_inc=sumID$ID_SITE[sumID$NUM_SAMP>=20]
ID=pfSiteSel(ID=sites_inc)

## Associated plots
plot(ID,zoom="sites")

## Simple test for transforming data
# Select site 1 (Cygnet Lake)

ID=pfSiteSel(ID=1)
plot(ID)

# Transformation to data
TR=pfTransform(ID,method=c("MinMax", "Box-Cox", "Z-Score"))

# Plot Transformed and raw data
# First retrieve raw data for Cygnet using pfExtract

RAW=pfExtract(ID=1)

dev.off()
par(mfrow=(c(2,1)))

plot(RAW[,3],RAW[,4],type="l")
plot(TR$Age,TR$TransData,type="l")


## Transforming and Compositing
## Example 1: Usage as in Power et al. 2008
## Data transformation
ID=pfSiteSel(Latlim=c(30,70),Longlim=c(-100,-50))
TR1=pfTransform(ID, method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,2000))

## Diagnostic pdf file with transformed series:
pfDiagnostic(ID, method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,2000))


## Compositing: basic binning procedure
COMP=pfComposite(TR1, binning=TRUE, bins=seq(0,12000,500))
plot(COMP)

## The resulting matrix can be saved
write.csv(COMP$Result,file="temp.csv")


## Compositing: Using the locfit package equivalent procedure to Daniau et al. 2012

COMP2=pfCompositeLF(TR1, tarAge=seq(-50,12000,20), binhw=20, hw=500,nboot=100)
plot(COMP2)

## And save
write.csv(COMP2$Result,file="temp2.csv")


## Example 3; Circular block bootstrapp
COMP=pfComposite(TR1, binning=TRUE, bins=seq(0,2000,100))
sea1=pfCircular(COMP,b=3,conf=c(0.005,0.025,0.975,0.995),nboot=100)

# Figure 
plot(sea1)
}
