
pfToKml=function(x,file="NULL")
{
  if(file=="NULL") stop("Output not specified")

## For testing
# x=pfSiteSel(lat> 40, long>-85, long<(-70), lat<60)
# plot(x) ...
# file="/Users/Olivier/Desktop/truc.kml"

df=data.frame(name=row.names(summary(x)),summary(x))
sp::coordinates(df) = ~long+lat

sp::proj4string(df)<-sp::CRS("+init=epsg:4326")

options(warn=-1)
rgdal::writeOGR(df, dsn=file, layer= "df", driver="KML", dataset_options=c("NameField=name"))

}