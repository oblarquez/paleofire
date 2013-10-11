install.packages("devtools")
library(devtools)

install_github("paleofiRe","oblarquez")
library(paleofire)

## Simple test

ID=pfSiteSel(Region="ASIA")
plot(ID)

