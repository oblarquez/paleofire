paleofire
=========
paleofire: an R package to analyse sedimentary charcoal records from the Global Charcoal Database to reconstruct past biomass burning

The `paleofire` package provides tools to extract and analyse charcoal sedimentary data stored in the Global Charcoal Database (see http://gpwg.paleofire.org for details). Main functionalities includes data extraction and sites selection, transformation and homogenization of the charcoal records as well as regional to global compositing.


Installation:
=============

To install `paleofire` from CRAN simply type `install.packages("paleofire")` at the R prompt.  

To install `paleofire` from GitHub the `devtools` package is required: on Windows platform the Rtools.exe program is required in order to install the `devtools` package. Rtools.exe can be downloaded for a specific R version on http://cran.r-project.org/bin/windows/Rtools/

Once `devtools` is installed type the following lines at R prompt to install `paleofire`:

```R
library(devtools)
install_github(repo="paleofire", username="paleofire", ref="master")
# The GCD package is required along paleofire:
install_github(repo="GCD", username="paleofire", ref="master")
library(paleofire)
```

To test everything is working you can plot a map of all charcoal records included in the Global Charcoal Database:

```R
plot(pfSiteSel())
```

For details and examples about `paleofire` please refer to the included [manual](http://cran.r-project.org/web/packages/paleofire/paleofire.pdf) or [vignette](http://cran.r-project.org/web/packages/paleofire/vignettes/paleofire-paper.pdf).

Citation:
=============

To cite `paleofire` in publications use:

  Blarquez O., Vannière B., Marlon J. R., Daniau A.-L. Power M. J. Brewer S.
  and Bartlein P. J. 2014. paleofire: An R package to analyse sedimentary
  charcoal records from the Global Charcoal Database to reconstruct past
  biomass burning. Computers and Geosciences, 72, 255-261
  
```R
citation (package="paleofire")
```

Maintainer: Olivier Blarquez <blarquez@gmail.com>
