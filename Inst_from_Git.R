install.packages("devtools")
library(devtools)

install_github("paleofire","paleofire")
library(paleofire)

## Simple test

ID=pfSiteSel(Region="ASIA")
plot(ID)

