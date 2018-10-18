#' Superposed Epoch Analysis
#' 
#' The function performs a Superposed Epoch Analysis on a provided temporal serie. 
#' The function uses pfCircular function for the computation of 
#' the block bootstrap procedure. The function could be used on both dendrochronological 
#' data and other data expressed in CE ages as well on paleoecological series expressed in BP.
#' Proxy data and ages must be spaced evenly but not necessarily using 1 yr time steps. 
#'
#' @param x data frame or matrix with ages and proxy values, younger ages on top. 
#' @param y events dates.
#' @param lag lag time used for calculating the SEA.
#' @param b A numeric giving block size, if NULL the optimal block size is given by: b= 2x(-1 /log(p)), where p is the lag one
#' autocorrelation coefficient of the serie (Adams, Mann & Ammann 2003).
#' @param conf confidence intervals for the block bootstrap procedure.
#' @param nboot number of bootstrap replicates.
#' @param age type of ages used in x[,1] either "CE" for Common Era or "BP" for Before Present.
#'
#' @return
#' \item{res}{A "pfCircular" object with estimated confidence intervals.}
#' 
#' @export
#' @author O. Blarquez
#' @seealso \code{\link{pfCircular}}
#' @examples
#' \dontrun{
#' ## Generate some fake data
#' set.seed(1)
#' n <- 100 # number of data points
#' t <- seq(0,4*pi,,100)
#' a <- 3
#' b <- 2
#' c.unif <- runif(n)
#' amp <- 4
#'
#' # generate data and calculate "y"
#' set.seed(1)
#' y1 <- a*sin(b*t)+c.unif*amp # add uniform error
#'
#' # SEA applied to fake dendrochronological data in CE
#' plot(rev(seq(1901,2000,1)), y1, t="l", ylim=range(y1)*c(1.2))
#' y=c(1923,1948,1972,1995)
#' points(y,rep(0,length(y)))
#' x=data.frame(rev(seq(1901,2000,1)),value=y1)
#' lag=10
#'
#' #Perform SEA
#' res=SEA(x, y, lag, b = NULL, conf = c(0.05, 0.95), nboot = 1000, age="CE")
#' plot(res,xlim=c(-10,10),xlab="lag",ylab="Composite mean")
#' 
#' # SEA applied to fake paleoecological data in BP
#' plot(seq(-50,49,1), y1, t="l", ylim=range(y1)*c(1.2),xlim=c(50,-50))
#' y=1950-c(1923,1948,1972,1995)
#' points(y,rep(0,length(y)))
#' x=data.frame(seq(-50,49,1),value=y1)
#' # Perform SEA
#' res=SEA(x, y, lag, b = NULL, conf = c(0.05, 0.95), nboot = 1000, age="BP")
#' plot(res,xlim=c(-10,10),xlab="lag",ylab="Composite mean")
#'}


SEA=function(x, y, lag, b = NULL, conf = c(0.05, 0.95), nboot = 1000, age="CE"){
  
  x=as.matrix(x)
  
  if (is.null(b) == TRUE) {
    r <- cor(x[1:length(x[, 1]) - 1, 2], x[2:length(x[, 1]), 2], use = "pairwise.complete.obs")
    b <- ceiling(2 * (-1 / log(abs(r))))
    b[b == 0 | b == 1 | is.na(b) | is.finite(b) == FALSE] <- 2
  }
  
  tunit=abs(mean(diff(x[,1])))
  x1=cbind(seq(from=min(x[,1])-lag*tunit , to=min(x[,1])-tunit, by=tunit),rep(NA,lag))
  x2=cbind(seq(from=max(x[,1])+tunit , to=max(x[,1])+lag*tunit, by=tunit),rep(NA,lag))
  
  if(all(x[,1] == cummax(x[,1])) & age=="BP"){
    x=rbind(x1,x,x2) 
  } 
  if(all(x[,1] == cummax(x[,1]))==FALSE & age=="CE"){
    x=rbind(apply(x2, 2, rev), x,apply(x1, 2, rev))
  } 
  if(all(x[,1] == cummax(x[,1])) & age=="CE"){
    stop("wrong age order")
  } 
  if(all(x[,1] == cummax(x[,1]))==FALSE & age=="BP"){
    stop("wrong age order")
  } 
  
  
  ou=which(x[,1] %in% y)  
  M=matrix(ncol=length(y), nrow=lag*2+1)  
  
  for(i in 1:length(y)){
    M[,i]=x[seq(from=ou[i]-lag, to=ou[i]+lag, by=1),2]
  }
  
  obj=list()
  obj$BinnedData=M
  obj$BinCentres=seq(-lag,+lag,1)
  res=pfCircular(obj, b = b, conf = conf, nboot = nboot)
  res$BootMean=rowMeans(M,na.rm = T)
  res$yr=rev(seq(-lag,+lag,1))
  class(res)="pfCircular"
  return(res)
  
}



