# runShinyFMBN v 4_1 ------------------------------------------------------

# This script will:
# a. check if you have all the packages needed to run ShinyFMBN, 
#    and install and load them as needed
# b. run ShinyFMBN v 3
# c. set the working directory to the location of the script

# install/load packages ---------------------------------------------------

.cran_packages <-
  c(
    "shiny",
    "data.table",
    "DT",
    "tidyverse",
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
# I am using pak because it handles dependencies better
if(!require(pak, quietly = T)) {
  install.packages("pak")
  require(pak)
  }
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)) {
  pak::pkg_install(.cran_packages[!.inst])
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  pak::pkg_install(.bioc_packages[!.inst])
}
sapply(c(.cran_packages, .bioc_packages), require,
       character.only = TRUE)


runApp("ShinyFMBN_4")



# Acknowledgements

# This work was was carried out within the PRIN 2022 PNRR Project 
# NCY diversity P20229JMMH and received funding from the European Union 
# Next-GenerationEU, CUP C53D23007560001 (PIANO NAZIONALE DI RIPRESA E 
# RESILIENZA (PNRR) – MISSIONE 4 COMPONENTE 2,  INVESTIMENTO 1.4 – D.D. 
# 1364/2023). This script and its contents reflects only the authors’ 
# views and opinions,  neither the European Union nor the European 
# Commission can be considered  responsible for them.

# Copyright and license

# Copyright 2024, 2025, 2026 Eugenio Parente
# Permission is hereby granted, free of charge, to any person obtaining 
# a copy of this software and associated documentation files (the "Software"), 
# to deal in the Software without restriction, including without limitation 
# the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the Software 
# is furnished to do so, subject to the following conditions:
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
# SOFTWARE.


