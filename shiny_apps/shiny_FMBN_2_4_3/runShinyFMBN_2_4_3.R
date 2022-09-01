# runShinyFMBN v 2_4_3 ------------------------------------------------------

# This script will:
# a. check if you have all the packages needed to run ShinyFMBN, 
#    and install and load them as needed
# b. run ShinyFMBN v 2


# install/load packages ---------------------------------------------------

.cran_packages <-
  c(
    "shiny",
    "DT",
    "tidyverse",
    "reshape2",
    "stringr",
    "lubridate",
    "igraph",
    "magrittr",
    "randomcoloR",
    "forcats"
  )
.bioc_packages <- c(
  "BiocManager", 
  "phyloseq"
)
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  if(!.inst[1]) {
    install.packages("BiocManager")
    .inst <- .bioc_packages %in% installed.packages()
  }
  if(any(!.inst[2:length(.inst)])) {
    BiocManager::install(.bioc_packages[!.inst], ask = F)
  }
}
sapply(c(.cran_packages, .bioc_packages), require,
       character.only = TRUE)


runApp("ShinyFMBN2_4_3")

# not run
# in some cases installation of Bioconductor packages may fail
# try this
# installing core bioconductor packages -----------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# check R version
if (as.numeric(R.version$major) >= 4) {
  BiocManager::install(version = "3.15")
} else {
  BiocManager::install(version = "3.10")
}
