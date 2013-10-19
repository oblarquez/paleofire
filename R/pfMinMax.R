pfMinMax=function(serie){
  # Minimax transformation
serie <- (serie - min(serie,na.rm=TRUE))/(max(serie,na.rm=TRUE)-min(serie,na.rm=TRUE))
}