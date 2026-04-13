# run ShinyFMBN from GitHub
if(!require(shiny)) install.packages("shiny")
require(shiny)
shiny::runGitHub(
  repo = "FoodMicrobionet",
  username = "ep142",
  subdir = "shiny_apps/shiny_FMBN_4_light/ShinyFMBN_4"
)
