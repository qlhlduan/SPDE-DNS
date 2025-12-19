## ===============================
## Package requirements
## ===============================

required_packages <- c(
  "INLA",
  "rSPDE",
  "fmesher",
  "inlabru",
  "Matrix",
  "sf",
  "parallel",
  "fields",
  "reshape2",
  "ggplot2",
  "gstat",
  "quadprog",
  "dplyr",
  "reshape2"
)

## Install missing packages
installed <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

## Load packages
invisible(lapply(required_packages, library, character.only = TRUE))

## INLA-specific repository (important!)
if (!"INLA" %in% installed) {
  install.packages(
    "INLA",
    repos = c(getOption("repos"),
              INLA = "https://inla.r-inla-download.org/R/stable"),
    dep = TRUE
  )
}
