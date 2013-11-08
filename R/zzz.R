.onAttach <- function(lib, pkg) {
  packageStartupMessage("This is paleofire v",utils::packageDescription("paleofire",field="Version"),appendLF = TRUE)
  checkGCDversion()
} 

checkGCDversion <- function() {
  # Check to see if installed
  if (!"GCD" %in% utils::installed.packages()[, 1]) {
    Checks <- "Failed"
  } else {
    # Compare version numbers
    temp <- getURL("https://raw.github.com/paleofire/GCD/daily/DESCRIPTION")      
    CurrentVersion <- gsub("^\\s|\\s$", "", 
                           gsub(".*Version:(.*)\\nDate.*", "\\1", temp))
    
    if (utils::packageVersion("GCD") == CurrentVersion) {
      Checks <- "Passed"
    }
    if (utils::packageVersion("GCD") < CurrentVersion) {
      Checks <- "Failed"
    }
  }
  
  switch(
    Checks,
    Passed = { message("Everything looks OK! GCD up to date: v",CurrentVersion) },
    Failed = {
      ans = readline(
        "'GCD is either outdated or not installed. Update now? (y/n) ")
      if (ans != "y")
        return(invisible())
      install_github("GCD",username="paleofire",ref="daily")
    })
  # Some cool things you want to do after you are sure the data is there
  library(GCD)
}


## Test for daily branch daily


