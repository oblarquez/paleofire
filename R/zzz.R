.onAttach <- function(lib, pkg) {
  #checkGCDversion()
  packageStartupMessage("This is paleofire v",utils::packageDescription("paleofire",field="Version"),appendLF = TRUE)
} 

checkGCDversion <- function() {
  # Check to see if installed
  if (!"GCD" %in% utils::installed.packages()[, 1]) {
    Checks <- "Failed"
  } else {
    # Compare version numbers
    temp <- getURL("https://raw.github.com/paleofire/GCD/master/DESCRIPTION",
                   ssl.verifypeer = FALSE)      
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
#       ans = readline(
#         "GCD is either outdated or not installed. Update now? (y/n) ")
#       if (ans != "y")
#         return(invisible())
      packageStartupMessage("GCD is either outdated or not installed. Installing...")
      install_github("GCD",username="paleofire",ref="master")
    })
  # packageStartupMessage("This is paleofire v",utils::packageDescription("paleofire",field="Version"),appendLF = TRUE)
}
