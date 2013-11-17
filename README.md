paleofire
=========
paleofire: an R package to analyse sedimentary charcoal records from the Global Charcoal Database to reconstruct past biomass burning

The `paleofire` package provides tools to extract and analyse charcoal sedimentary data stored in the Global Charcoal Database (see http://gpwg.org for details). Main functionalities includes data extraction and sites selection, transformation and homogenization of the charcoal records as well as regional to global compositing.


Installation:
=============
To install `paleofire` from the present GitHub repository the `devtools` package is required: on Windows platform the Rtools.exe program is required in order to install the `devtools` package. Rtools.exe can be downloaded for a specific R version on http://cran.r-project.org/bin/windows/Rtools/

Once `devtools` is installed type the following lines at R prompt to install `paleofire`:

```R
library(devtools)
install_github("paleofire","paleofire")
library(paleofire)
```

To test everything is working you can plot a map of all charcoal records included in the Global Charcoal Database:

```R
plot(pfSiteSel())
```

For details and examples about `paleofire` please refer to the included [manual](https://github.com/paleofire/paleofire/raw/master/paleofire-manual.pdf).

Maintainer: Olivier Blarquez <blarquez@gmail.com>
