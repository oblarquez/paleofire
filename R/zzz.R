.onAttach <- function(lib, pkg) {
  #checkGCDversion()
  packageStartupMessage("This is paleofire v",utils::packageDescription("paleofire",field="Version"),appendLF = TRUE)
} 

checkGCDversion <- function() {
  # Check to see if installed
  if (!"GCD" %in% utils::installed.packages()[, 1]) {
    Checks <- "Failed"
  } else {
    oldones=old.packages()
    if ("GCD" %in% oldones[,1]) {
      Checks <- "Failed"
    } else Checks <- "Passed"
  }
  switch(
    Checks,
    Passed = { message("Everything looks OK! GCD up to date: v", packageVersion("GCD")) },
    Failed = {
#       ans = readline(
#         "GCD is either outdated or not installed. Update now? (y/n) ")
#       if (ans != "y")
#         return(invisible())
      packageStartupMessage("GCD is either outdated or not installed. Installing...")
      install.packages("GCD")
    })
  # packageStartupMessage("This is paleofire v",utils::packageDescription("paleofire",field="Version"),appendLF = TRUE)
}
