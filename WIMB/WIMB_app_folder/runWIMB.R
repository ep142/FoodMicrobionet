# runWIMB ------------------------------------------------------

# This script will:
# a. check if you have all the packages needed to run the wimb_app, 
#    and install and load them as needed
# b. run wimb_app


# install/load packages ---------------------------------------------------

.cran_packages <-
  c(
    "shiny",
    "bslib",
    "bsicons",
    "markdown",
    "tidyverse",
    "data.table",
    "DT",
    "scales"
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
.loaded <- sapply(.cran_packages, require, character.only = TRUE)
if (any(!.loaded)) message("Failed to load: ", paste(names(.loaded)[!.loaded], collapse = ", "))

runApp("wimb_app", launch.browser = TRUE)


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


