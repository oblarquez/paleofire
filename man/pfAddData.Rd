\name{pfAddData}
\alias{pfAddData}
\title{
Add user defined charcoal data series to paleofire}
\description{
This function is used to create a "pfAddData" object, from user defined csv files containing charcoal data, to be passed to pfTransform. Usually csv files should contain three columns with Depth, Age, Charcoal quantity in this same order. A metadata csv file should also be specified with sites location information (three columns with: SITE_NAME, LATITUDE, LONGITUDE). CharAnalysis data files could also be used, in this case the file must include the following informations: DepthTop, DepthBottom, AgeTop, AgeBottom, Volume and Charcoal value in this exact order. Then the files are passed to the pretreatment function in order to calculate Charcoal Accumulation Rates (see pretreatment for details).    
}
\usage{
pfAddData(files, metadata, type = "NULL", Int = TRUE, first =
                 NULL, last = NULL, yrInterp = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{
Character, names and path to csv files.}
 \item{metadata}{
Character, names and path to metadata csv file.}
 \item{type}{
Character, "NONE": user defined csv (default), "CharAnalysis": CharAnalysis data file.}
  \item{Int}{
Logical specifying whether the pretreatment function interpolates particle zero counts, default TRUE.}
  \item{first,last}{
Numeric, date of the first, last sample for accumulation rate calculation, if NULL first, last are automatically specified as the the minimum and maximum ages of the record respectively.}
  \item{yrInterp}{
Numeric, temporal resolution of the interpolated accumulation rates, if NULL, yrInterp is automatically specified as the median resolution of the record.}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{out}{A list with merged data files that can be passed to \code{\link{pfTransform}}}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
O. Blarquez}

\seealso{\code{\link{pretreatment}}}

\examples{
# not run

## Ad user own data from CharAnalysis file (csv)
## In this example we will use data from:
# Senici, D., A. Lucas, H. Y. H. Chen, Y. Bergeron, A. Larouche, B. Brossier, O.
#Blarquez, and A. A. Ali. 2013. Multi-millennial fire frequency and tree abundance
#differ between xeric and mesic boreal forests in central Canada. Journal of Ecology:
#101, 356-367.

# files=c("http://blarquez.com/public/data//Ben.csv",
#        "http://blarquez.com/public/data/Small.csv")
# metadata=c("http://blarquez.com/public/data/metadata.csv")

# mydata=pfAddData(files=files,metadata=metadata,type="CharAnalysis")

## Transform and compositing:
# TR1=pfTransform(add=mydata, method=c("MinMax","Box-Cox","Z-Score"),
#                 BasePeriod=c(200,2000))
# COMP2=pfCompositeLF(TR1, tarAge=seq(-50,8000,20), hw=500, nboot=100)
# plot(COMP2)
}