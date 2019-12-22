
# runShinyFMBN ------------------------------------------------------------

# This script will:
# a. check if you have all the packages needed to run ShinyFMBN, 
#    and install and load them as needed
# b. run ShinyFMBN v 1_2

# install/load packages -------------------------------------------------

.cran_packages <- c("shiny", "DT", "tidyverse", "reshape2", "stringr", 
                    "tibble", "lubridate" ,"igraph")
.bioc_packages <- c("BiocManager", "phyloseq")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  if(!"BiocManager" %in% installed.packages()) install.packages("BiocManager")
  BiocManager::install(.bioc_packages[2:length(.inst)])
}
sapply(c(.cran_packages, .bioc_packages), require, 
       character.only = TRUE) 

runApp("ShinyFMBN1_2")
