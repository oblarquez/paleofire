.onLoad <- function(lib, pkg) {
  packageStartupMessage("This is paleofire v",utils::packageDescription("paleofire",field="Version"),appendLF = TRUE ," (GCDv03)")
} 