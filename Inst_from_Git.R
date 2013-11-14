install.packages("devtools")
library(devtools)
remove.packages("GCD")
remove.packages("paleofire")

install_github("paleofire","paleofire")
library(paleofire)

## Simple test

ID=pfSiteSel(Region="ASIA")
plot(ID)

