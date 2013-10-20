CheckVersionFirst <- function() {
  # Check to see if installed
  if (!"GCD" %in% installed.packages()[, 1]) {
    Checks <- "Failed"
  } else {
    # Compare version numbers
    require(RCurl)
    temp <- getURL("https://raw.github.com/paleofire/GCD/master/DESCRIPTION")    
    CurrentVersion <- gsub("^\\s|\\s$", "", 
                           gsub(".*Version:(.*)\\nDate.*", "\\1", temp))
    
    if (packageVersion("GCD") == CurrentVersion) {
      Checks <- "Passed"
    }
    if (packageVersion("GCD") < CurrentVersion) {
      Checks <- "Failed"
    }
  }
  
  switch(
    Checks,
    Passed = { message("Everything looks OK! GCD data up to date") },
    Failed = {
      ans = readline(
        "'GCD is either outdated or not installed. Update now? (y/n) ")
      if (ans != "y")
        return(invisible())
      require(devtools)
      install_github("paleofire", "GCD")
    })
  # Some cool things you want to do after you are sure the data is there
}