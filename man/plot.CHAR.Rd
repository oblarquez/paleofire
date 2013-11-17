\name{plot.CHAR}
\alias{plot.CHAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot CHAR}
\description{
Plot an object of the class "CHAR" returned by the pretreatment function.
Original accumulation rates are presented using grey bars, accumulation rates interpolated at equal time steps are presented by a black curve.}

\usage{
\method{plot}{CHAR}(x,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of the class "CHAR".}
\item{\dots}{\dots}
}

\author{
O. Blarquez}

\examples{
## In this example we will use the charcoal record of the Lac du Loup (Blarquez et al. 2010)
## Load raw charcoal data in mm^2
A=read.csv("http://blarquez.com/public/code/loupchar.csv")
C_=A[,6] # charcoal areas
P_=A[,1:5] # CmTop, CmBot, AgeTop, AgeBot, Volume 


## Calculates charcoal accumulation rate (CHAR, mm2.cm-2.yr-1)
CHAR=pretreatment(params=P_,serie=C_,Int=TRUE)
plot(CHAR)

}
