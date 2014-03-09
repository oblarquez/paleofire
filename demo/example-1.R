## Interactive sites selection:
# ID=pfInteractive()
library(devtools)
install_github("paleofire", repo="GCD",ref="release3.0.1")
library(GCD)
library(paleofire)
## Site selection using criterions
# DateInt parameter is used to set the mean interval which is required between two
# dating points (ex 14C) for sites to be selected for a complete list of criterions
# that can be used see pfSiteSel function


ID=pfSiteSel(lat>40)

plot(ID)

ID1=pfSiteSelM(ID_SITE %in% ID$SitesIDS & LONGITUDE>-90 & LONGITUDE<(-50))
plot(ID1,xlim=c(-100,-50),ylim=c(30,70))

plot(ID1,zoom="world")


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

# Transformation of data
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
ID=pfSiteSel(Latlim=c(50,70),Longlim=c(-90,-50), DateInt=3000)
TR1=pfTransform(ID, method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,2000))

## Diagnostic pdf file with transformed series:
# pfDiagnostic(ID, method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,2000), 
# FileName = "Diagnostic.pdf")


## Compositing: basic binning procedure
COMP=pfComposite(TR1, binning=TRUE, bins=seq(0,8000,500))
plot(COMP)

## The result matrix can be saved
write.csv(COMP$Result,file="temp.csv")


## Compositing: Using the locfit package equivalent procedure to Daniau et al. 2012

COMP2=pfCompositeLF(TR1, tarAge=seq(-50,8000,20), binhw=20, hw=500,nboot=100)
plot(COMP2)

## And save