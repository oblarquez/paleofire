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

ID1=pfSiteSel(id_site==1)
plot(ID1)

# Transformation of data
TR=pfTransform(ID1,method=c("MinMax", "Box-Cox", "Z-Score"))

# Plot Transformed and raw data
# First retrieve raw data for Cygnet using pfExtract

RAW=pfExtract(ID1)

dev.off()
par(mfrow=(c(2,1)))

plot(RAW[,3],RAW[,4],type="l")
plot(TR$Age,TR$TransData,type="l")


## Transforming and Compositing
## Example 1: Usage as in Power et al. 2008
## Data transformation
TR1=pfTransform(ID, method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,2000),MethodType="SIEV")

TR1=pfTransform(ID, method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,2000))

## Diagnostic pdf file with transformed series:
# pfDiagnostic(ID, method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,2000), 
# FileName = "Diagnostic.pdf")


## Compositing: basic binning procedure
COMP=pfComposite(TR1, binning=TRUE, bins=seq(0,8000,500))
plot(COMP)

## The result matrix can be saved
# write.csv(COMP$Result,file="temp.csv")


## Compositing: Using the locfit package equivalent procedure to Daniau et al. 2012

COMP2=pfCompositeLF(TR1, tarAge=seq(-50,8000,20), binhw=20, hw=500,nboot=100)
plot(COMP2)



## Composite charcoal record for western North America:
ID=pfSiteSel(id_region=="WNA0" & l12==1)

plot(ID, zoom="world" )
## Transform data
res3=pfTransform(ID,method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,4000),Interpolate=FALSE)

## Composite
comp=pfComposite(res3,bins=seq(from=-500,to=12500,by=1000))
plot(comp)

## Kruskal Wallis Anova
comparison=pfKruskal(comp,alpha=0.05, p.adj="none")
plot(comparison,trend=TRUE)

## Kruskal Wallis Anova
comparison=pfKruskal(res3,alpha=0.05,bins=seq(from=-500,to=12500,by=1000))
plot(comparison,ylim=c(-5,7))


## Here

ID=pfSiteSel(id_region=="ENA0", l12==1, long>-85)

res3=pfTransform(ID,method=c("MinMax","Box-Cox","Z-Score"),BasePeriod=c(200,4000))

p=pfGridding(res3,age=1000)
summary(p)

plot(p,points=T,empty_space=100)




