.onLoad <- function(libname, pkgname) {
  paramDir <- system.file("PAM50/bioclassifier_R", package = pkgname)
  source(file.path(paramDir, "subtypePrediction_functions_PNmodified.R"))
  source(file.path(paramDir, "subtypePrediction_distributed_PNmodified.R"))
}
