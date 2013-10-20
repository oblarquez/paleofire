.onAttach <- function(lib, pkg) {
  packageStartupMessage("This is paleofire v",utils::packageDescription("paleofire",field="Version"),appendLF = TRUE)
} 
.onLoad <- function(lib, pkg) {
  CheckVersionFirst()
} 