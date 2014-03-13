\name{pfSiteSel}
\alias{pfSiteSel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GCD sites selection methods}
\description{
Main function used for site selection, uses data strored in data(paleofiresites) to perform site selection according to multiple criterion, those criterions could be either geographic, based on series attributes (e.g. # of datings), or on sites attributes (e.g. biome).}
\usage{
pfSiteSel(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{...}{Any combination of conditions defined by relational operators and or logical operators that are applied on the "paleofiresites" dataset. See examples below:}
}
    
\details{
Use data(paleofiresites);names(paleofiresites) to retrieve the conditions that could be used to select sites i.e.:
 id_site, site_name, lat, long, elev, pref_units, biome, id_region, id_country, id_site_type, water_depth, 
 id_basin_size, id_catch_size, id_land_desc, dating_type, min_est_age, max_est_age, num_dating, age_model, 
 data_source, qtype, rf99, l12, num_samp, date_int.}
\value{
An object of the class "pfSiteSel" (list) with "id_site" and "site_name" components.
}
\seealso{\code{\link[GCD]{paleofiresites}}}
\author{
O. Blarquez}

\examples{
## Sites selection examples

## Select all sites
ID=pfSiteSel()

## Site in the Biome #8
ID=pfSiteSel(biome==8)
plot(ID,zoom="world")

## Site in the biome #8 or in the biome #6
ID=pfSiteSel(biome==8 | biome==6)

## Sites in North America by geographic location
ID=pfSiteSel(lat>25, lat<75, long<(-45), long>-150) 
plot(ID,zoom="world")

## is equivalent to:
ID=pfSiteSel(lat>25 & lat<75 & long<(-45) & long>-150) 
plot(ID,zoom="world")

## By region criterion
ID=pfSiteSel(id_region==c("ENA0","WNA0"))
plot(ID,zoom="world")

## WRONG, use the %in% operator when concatenating two characters
# ID=pfSiteSel(id_region %in% c("ENA0","WNA0"))
# plot(ID,zoom="world")

## Pas-de-Fond site
pfSiteSel(site_name=="Pas-de-Fond")

## All sites in  eastern North America that are not Pas-de-Fond
pfSiteSel(site_name!="Pas-de-Fond", id_region=="ENA0")

## Sites with on average one dating point every 250 to 300 yrs
pfSiteSel(date_int>=250 & date_int<=300)

## Sites between 0, 100 m elevation in Asia
ID=pfSiteSel(elev>0 & elev<100, id_region=="ASIA")

## All sites that are not marine nor fluvial
ID=pfSiteSel(id_land_desc!="MARI" , id_site_type!="FLUV" & id_site_type!="LFLU")
plot(ID)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
