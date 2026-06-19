# run WIMB_app from GitHub v.0.1 19/06/26

options(timeout = 600)
if(!require(shiny)) install.packages("shiny")
require(shiny)
shiny::runGitHub(
  repo = "FoodMicrobionet",
  username = "ep142",
  subdir = "WIMB/WIMB_app_folder/WIMB_app",
  launch.browser = T
)
