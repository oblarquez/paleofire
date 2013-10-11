pfMinMax=function(serie){
  # Minimax transformation
serie <- (serie - min(serie))/(max(serie)-min(serie))
}