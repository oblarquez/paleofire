## Interactive sites selection:
# ID=pfInteractive()


library(devtools)
install_github("paleofire", repo="GCD",ref="release3.0.1")
install_github("paleofire", repo="paleofire",ref="release1.1")

library(GCD)
library(paleofire)

## Site selection using criterions

data(paleofiresites)
names(paleofiresites)

ID=pfSiteSel(id_land_desc!="MARI" , id_site_type!="FLUV" & id_site_type!="LFLU")
plot(ID)

ID=pfSiteSel(lat>0 & rf99==6)
plot(ID)

ID=pfSiteSel(lat>0, rf99==6 | l12==1)
plot(ID)

ID=pfSiteSel(lat>0, rf99==6 | l12==1)
plot(ID)

ID=pfSiteSel(lat>0, rf99==6 | l12==1, date_int<=2000 & num_samp>30)
plot(ID)

ID=pfSiteSel(lat>0, rf99==6 | l12==1, elev<=1000)
plot(ID)

## Filter sites based on sample number using summary function


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