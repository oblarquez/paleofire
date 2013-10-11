\name{paleofiresites}
\alias{paleofiresites}
\docType{data}
\title{
paleofiresites}
\description{
Sites metadata associated with paleofiredata dataset}
\usage{data(paleofiresites)}
\format{
  A data frame with 736 observations on the following 21 variables.
  \describe{
    \item{\code{ID_SITE}}{a numeric vector}
    \item{\code{SITE_NAME}}{a character vector}
    \item{\code{LATITUDE}}{a numeric vector}
    \item{\code{LONGITUDE}}{a numeric vector}
    \item{\code{ELEV}}{a numeric vector}
    \item{\code{PREF_UNITS}}{a factor with levels \code{125} \code{DWT} \code{TOF} \code{WTO} \code{100P} \code{125V} \code{5015} \code{7PRA} \code{8POS} \code{ABSC} \code{ARCC} \code{AREA} \code{ARPG} \code{BCCT} \code{C010} \code{C100} \code{C120} \code{C125} \code{C140} \code{C150} \code{C250} \code{C550} \code{CARE} \code{CCMY} \code{CG2M} \code{CHCO} \code{CHRX} \code{CM10} \code{CMGR} \code{CMML} \code{COML} \code{CP20} \code{CPAR} \code{CPRA} \code{DC13} \code{FRAG} \code{G05I} \code{GCHG} \code{GCMY} \code{GRGR} \code{GRPC} \code{GT05} \code{GT10} \code{GT12} \code{GT15} \code{GT18} \code{GT24} \code{GT25} \code{GT2X} \code{GT30} \code{GT38} \code{GT50} \code{GT53} \code{GT5C} \code{GT5D} \code{LT12} \code{LT15} \code{LT18} \code{LT1M} \code{LT25} \code{LT40} \code{M2GR} \code{M2ML} \code{MCM3} \code{MM2C} \code{MM2G} \code{MMDW} \code{MMML} \code{MT10} \code{ORDI} \code{OTHE} \code{PCM1} \code{PCMY} \code{PHYT} \code{PP25} \code{PP51} \code{PP52} \code{PPOL} \code{PRAB} \code{PROB} \code{PTC0} \code{PTC2} \code{PTC3} \code{PTCT} \code{SQC0} \code{SQCU} \code{SQCY} \code{SQG5} \code{TOM1} \code{TOM2} \code{TOM3} \code{TOM4} \code{TOM5} \code{TOM6} \code{TOM7} \code{TOM9} \code{TOMX} \code{V125} \code{X01K} \code{X105} \code{X106} \code{X120} \code{X125} \code{X150} \code{X160} \code{X180} \code{X18C} \code{X20P} \code{X250} \code{X25P} \code{X37P} \code{X459} \code{X512} \code{X520} \code{X550} \code{X55P} \code{XA1K} \code{XARE} \code{XARP} \code{XARS} \code{XC25} \code{XCM3} \code{XCOP} \code{XFML} \code{XFRG} \code{XFRP} \code{XFRS} \code{XIMG} \code{XPEA} \code{XPEE} \code{XPIX} \code{XRIF}}
    \item{\code{BIOME}}{a numeric vector}
    \item{\code{ID_REGION}}{a factor with levels \code{AFRI} \code{ASIA} \code{AUST} \code{CENA} \code{ENA0} \code{EURO} \code{RUSS} \code{SSA0} \code{TROP} \code{TSA0} \code{WNA0}}
    \item{\code{ID_COUNTRY}}{a factor with levels \code{ALGE} \code{ARGE} \code{AUST} \code{BOLI} \code{BRAZ} \code{BULG} \code{BURU} \code{CANA} \code{CHIL} \code{CHIN} \code{COLO} \code{COST} \code{CROA} \code{DENM} \code{DOMR} \code{ECUA} \code{EL S} \code{ETHI} \code{FIJI} \code{FINL} \code{FRAN} \code{GEOR} \code{GERM} \code{HAIT} \code{HOND} \code{HUNG} \code{INDI} \code{INDO} \code{ITAL} \code{JAPA} \code{JORD} \code{KENY} \code{KYRG} \code{MADA} \code{MEXI} \code{MONG} \code{NEWC} \code{NEWZ} \code{NICA} \code{NIGE} \code{NORW} \code{OCEA} \code{PALA} \code{PANA} \code{PAPU} \code{PERU} \code{RUSS} \code{SLOV} \code{SLVE} \code{SOAF} \code{SOLO} \code{SPAI} \code{SWED} \code{SWIT} \code{SYRI} \code{TANZ} \code{THAI} \code{TRKY} \code{UGAN} \code{UK00} \code{UKRA} \code{US00} \code{VANU}}
    \item{\code{ID_SITE_TYPE}}{a factor with levels \code{CEST} \code{COAS} \code{CSLT} \code{FLUV} \code{HOLL} \code{LACU} \code{LCIR} \code{LDRA} \code{LEXP} \code{LFLU} \code{LGLA} \code{LGMC} \code{LKET} \code{LMOR} \code{LNAT} \code{LPLY} \code{LSAL} \code{LTEC} \code{LTEK} \code{LTHK} \code{LVOL} \code{MARI} \code{NOTK} \code{OTHE} \code{TBLA} \code{TBOG} \code{TDET} \code{TERR} \code{TFBU} \code{TFEN} \code{TFLO} \code{TFLP} \code{TGSW} \code{TLIT} \code{TLTP} \code{TLUS} \code{TMIR} \code{TMOR} \code{TMSH} \code{TMSP}}
    \item{\code{WATER_DEPTH}}{a numeric vector}
    \item{\code{ID_BASIN_SIZE}}{a factor with levels \code{} \code{LARG} \code{MED} \code{MEDI} \code{NOTA} \code{NOTK} \code{SMAL} \code{VEGL} \code{VELG} \code{vESM} \code{VESM}}
    \item{\code{ID_CATCH_SIZE}}{a factor with levels \code{} \code{LARG} \code{MED} \code{MEDI} \code{NOTA} \code{NOTK} \code{Smal} \code{SMAL}}
    \item{\code{ID_LANDS_DESC}}{a factor with levels \code{} \code{BARA} \code{COAS} \code{COPL} \code{DIUP} \code{DUNE} \code{HUDI} \code{ISUP} \code{LOPL} \code{LORI} \code{LRVA} \code{MARI} \code{MEAD} \code{MOUN} \code{NOTK} \code{OTHE} \code{ROUP} \code{UPPL}}
    \item{\code{DATING_TYPE}}{a factor with levels \code{} \code{AMS0} \code{CALO} \code{CALS} \code{DCAL} \code{NOTK} \code{OTHE} \code{RADI}}
    \item{\code{MIN_EST_AGE}}{a numeric vector}
    \item{\code{MAX_EST_AGE}}{a numeric vector}
    \item{\code{NUM_DATING}}{a numeric vector}
    \item{\code{AGE_MODEL}}{a character vector}
    \item{\code{DATA_SOURCE}}{a character vector}
    \item{\code{QTYPE}}{a character vector}
  }
}

\source{
GCD v3}

\examples{
data(paleofiresites)
}
\keyword{datasets}
