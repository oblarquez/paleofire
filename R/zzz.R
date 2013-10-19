.onAttach <- function(lib, pkg) {
 packageStartupMessage("This is paleofire v",utils::packageDescription("paleofire",field="Version"),appendLF = TRUE)
} 
# .onLoad <- function(lib, pkg) {
# #    fileURLsites <- "http://blarquez.com/public/data/paleofiresites.rda"
# #    #setInternet2(TRUE)
# #    download.file(fileURLsites ,destfile="./data/paleofiresites.rda",method="auto")
# #    fileURLdata <- "http://blarquez.com/public/data/paleofiredata.rda"
# #    #setInternet2(TRUE)
# #    download.file(fileURLdata ,destfile="./data/paleofiredata.rda",method="auto")
# # 
# #    gcd <- "http://blarquez.com/public/data/GCDver.rda"
# #    #setInternet2(TRUE)
# #    download.file(gcd ,destfile="./data/GCDver.rda",method="auto")
# #    load("./data/GCDver.rda")
# #    packageStartupMessage("Downloaded and installed GCD v",GCDver$ver, ",", GCDver$date)
# #    rm(GCDver)
# } 