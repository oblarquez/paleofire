
\name{pretreatment}
\alias{pretreatment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate particules accumulation rates for sediment records}
\description{
This is the R version of the CharAnalysis CharPretreatment.m function
originally develloped by P. Higuera and available at https://sites.google.com/site/charanalysis
}
\usage{
pretreatment(params, serie, Int = TRUE, first = NULL, last = NULL,
                 yrInterp = NULL)}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{serie}{
A proxy record to be transformed in accumulation rates, could be particule counts, surfaces, volumes, etc.}
  \item{params}{
A matrix with the following colums: CmTop, CmBot, AgeTop, AgeBot, Volume, in the same order.}
  \item{Int}{
Logical specifying whether the function interpolates particle zero counts, default TRUE}
  \item{first,last}{
Date of the first, last sample for accumulation rate calculation, if NULL first, last are automatically specified as the the minimum and maximum ages of the record respectively}
  \item{yrInterp}{
Temporal resolution of the interpolated accumulation rates, if NULL, yrInterp is automatically specified as the median resolution of the record}
}
\value{
Return an output structure with the following:
\item{cmI}{interpolated depths}
\item{ybpI}{interpolated ages}
\item{accI}{accumulation rates}
}

\author{
O. Blarquez translated from P. Higuera CharPretreatment.m function}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
# not run

## In this example we will use the charcoal record of the Lac du Loup from Blarquez et al. (2010).
# Blarquez, O., C. Carcaillet, B. Mourier, L. Bremond, and O. Radakovitch. 2010. Trees in the 
# subalpine belt since 11 700 cal. BP: origin, expansion and alteration of the modern forest. 
# The Holocene 20:139-146.

## Load raw charcoal data in mm^2
# A=read.csv("http://blarquez.com/public/code/loupchar.csv")
# C_=A[,6] # charcoal areas
# P_=A[,1:5] # CmTop, CmBot, AgeTop, AgeBot, Volume 


## Calculates charcoal accumulation rate (CHAR, mm2.cm-2.yr-1)
# CHAR=pretreatment(params=P_,serie=C_,Int=TRUE)
# plot(CHAR)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

