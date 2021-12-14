# assemble_FMBN, a script for assembling .rds files for FMBN

# v 4 02/12/2021
# this script takes as an input individual tables for Foodmicrobionet
# Studies
# Primers
# Samples
# FoodEx2exp
# Taxa
# Edges
# The tables must be in .xlsx formt (no formulas) and be in the same 
# folder. Creating a project allows to automate loading of the files.
# Version info has to be provided (see line 51)
# Once the tables have been imported in R several diagnostic checks are performed
# to identify potential mismatched in the tables. If the diagnostics fail no
# list is created.
# The tables, together with info on version, are assembled in two lists
# FMBN.rds which is compatible with ShinyFMBN 2.3
# FMBN_plus.rds which includes all the data but is not compatible with
# ShinyFMBN 2.3

# install/load packages

.cran_packages <- c("readxl","tidyverse")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}

sapply(c(.cran_packages), require, 
       character.only = TRUE)

# read filenames --------------------------------------------------------
filenames <- list.files()

# check that files with proper names are available
studies_file_flag <- "Studies.xlsx" %in% filenames
primers_file_flag <- "Primers.xlsx" %in% filenames
samples_file_flag <- "Samples.xlsx" %in% filenames
FoodEx2_file_flag <- "FoodEx2exp.xlsx" %in% filenames
taxa_file_flag <- "Taxa.xlsx" %in% filenames
edges_file_flag <- "Edges.xlsx" %in% filenames
file_flags <- c(studies_file_flag, primers_file_flag, samples_file_flag,
                FoodEx2_file_flag, taxa_file_flag,edges_file_flag)
if(all(file_flags)){
  cat("\nAll .xlsx files available\n")
} else {
  cat("\nOne or more .xlsx files are missing\n")
}

# version -----------------------------------------------------------------
version_n <- "4.1.2"
last_updated <- "02/12/2021"
version_text <- paste("This is version ", version_n, 
                      " of FoodMicrobionet. Last updated ", last_updated, ".", sep = "")
saveRDS(version_text, file = "version.rds")


# samples -----------------------------------------------------------------

# use the following commands if you want to load the file manually
# samples_fn <- file.choose()
sample_fn <- "Samples.xlsx"

samples <- read_xlsx(sample_fn,
                     col_types = c(
                       "text",
                       "numeric",
                       rep("text", 5),
                       rep("numeric", 3),
                       rep("text", 16)
                     ))

# check for duplicates

flag_dupli_sample <- !(sum(duplicated(samples$label_2))==0)

# find duplicated labels 
if(flag_dupli_sample){
  cat("WARNING: there are duplicated samples!!!!")
  dupli_sample_labels <- samples$label_2[duplicated(samples$label_2)]
}

# check if sample Id match row numbers
# this is necessary for sample extraction in the Shiny app: mismatches would 
# cause wrong samples to be selected in the aggregate tab

flag_sample_mismatch <- !(all(samples$sampleId == 1:nrow(samples)))

if(flag_sample_mismatch) {
  cat("WARNING: sampleIds DO NOT match row numbers")
  which(samples$sampleId != 1:nrow(samples))
  }

saveRDS(samples, file = "samples.rds")

# studies -----------------------------------------------------------------

# use the following command if you want to choose the file manually
# studies_fn <- file.choose()

studies_fn <- "Studies.xlsx"
studies <-
  read_xlsx(studies_fn,
            col_types = c(
              "numeric",
              rep("text", 5),
              "numeric",
              rep("text", 8),
              "numeric",
              rep("text", 12),
              rep("logical", 2)
            ))

# check if studyId in samples match those in studies
# set a flag
flag_study_mismatch <- !all(sort(unique(samples$studyId)) == sort(studies$studyId))
if(flag_study_mismatch){
  if(length(unique(samples$studyId))>nrow(studies)) {
    cat("WARNING: there are more studies in the sample table than in the study table")
    mismatched_studies <- setdiff(unique(samples$studyId), studies$studyId)
    cat("Mismatched studies: ", mismatched_studies)
  } else {
    cat("WARNING: there are more studies in the study table than in the sample table")
    mismatched_studies <- setdiff(studies$studyId, unique(samples$studyId))
    cat("Mismatched studies: ", mismatched_studies)
  }
}

saveRDS(studies, file = "studies.rds")


# primers -----------------------------------------------------------------

# use the following command if you want to choose the file manually
# primers_fn <- file.choose()

primers_fn <- "Primers.xlsx"
primers <- read_xlsx(primers_fn,
                     col_types = c(rep("text", 3),
                                   rep("numeric", 2),
                                   "text"))
saveRDS(primers, file = "primers.rds")

# FoodEx2 -----------------------------------------------------------------

# use the following command if you want to choose the file manually
# primers_fn <- file.choose()

foodex2_fn <- "FoodEx2exp.xlsx"
foodex2exp <- read_xlsx(foodex2_fn)
saveRDS(foodex2exp, file = "foodex2exp.rds")

# edges -------------------------------------------------------------------

# use the following command if you want to choose the file manually
# edges_fn <- file.choose()
edges_fn <- "Edges.xlsx"
# this is the largest file, import might take a while
edges <- read_xlsx(edges_fn)

# check if there are NA values and set a flag
flag_edges_na <- anyNA(edges)
if(flag_edges_na) {
  cat("WARNING: there are edges with NA values","\n")
  cat("The edges with NA values are at rows ", which(is.na(edges)))
  }

# check if sampleId in edges have the same length as those in samples

# set a flag when number of samples in samples and edges do not match
flag_edges_mismatch_1 <- !(length(unique(edges$sampleId)) == nrow(samples))
# initialize flag for mismatch_2
flag_edges_mismatch_2 <- ifelse(flag_edges_mismatch_1, T, F)

# check what is wrong
if(flag_edges_mismatch_1){
  cat("WARNING: the number of samples in the edges and in the sample table do not match","\n")
  if(length(unique(edges$sampleId))> nrow(samples)){
    cat("the excess edges are in the edge table at rows: ", setdiff(unique(edges$sampleId), samples$sampleId))
  } else {
    cat("the excess edges are in the sample table at rows: ", setdiff(samples$sampleId, unique(edges$sampleId)))
  }
} else {
  # check that sampleId are identical in both tables
  all(sort(unique(edges$sampleId)) == sort(samples$sampleId))
  flag_edges_mismatch_2 <- !all(sort(unique(edges$sampleId)) == sort(samples$sampleId))
  if(flag_edges_mismatch_2) cat("the mismatching edges are ", setdiff(unique(edges$sampleId), samples$sampleId))
}

saveRDS(edges, file = "edges.rds")


# taxa --------------------------------------------------------------------

# use the following command if you want to choose the file manually
# taxa_fn <- file.choose()
taxa_fn <- "Taxa.xlsx"
taxa <- read_xlsx(taxa_fn)
# check for duplicated names (should have no effect really)

# set a flag
flag_dupli_taxa <- !(sum(duplicated(taxa$label)) == 0)
# which are duplicated?
if(flag_dupli_taxa){
  dupli_taxa_labels <- taxa$label[duplicated(taxa$label)]
  cat("WARNING: there are duplicated taxa labels.","\n")
  cat("The duplicated taxa are ", dupli_taxa_labels)
}


# check that taxa match in the edges and taxa table
flag_taxa_mismatch <- !(all(unique(edges$taxonId) %in%  sort(taxa$taxonId)))
if(flag_taxa_mismatch){
  cat("WARNING: there are taxa in the edge table which are not found in the taxa table", "\n")
  missing_taxa <- which(!(unique(edges$taxonId) %in%  sort(taxa$taxonId)))
  cat("The taxonId for the missing taxa are: ", unique(edges$taxonId)[missing_taxa])
}


saveRDS(taxa, file = "taxa.rds")

references <- tibble(studyId = rep("Foodmicrobionet",3),
                        ref_complete = c(
                          paste0("Parente, E., Cocolin, L., De Filippis, F., Zotta, T., ",
                                 "Ferrocino, I., O’Sullivan, O., Neviani, E., De Angelis, ",
                                 "M., Cotter, P. D., Ercolini, D. 2016. FoodMicrobionet: a ",
                                 "database for the visualisation and exploration of food ",
                                 "bacterial communities based on network analysis. Int. J. ",
                                 "Food Microbiol. 219: 28-37.", sep =""),
                          paste0("De Filippis, F., Parente, E., Zotta, T., Ercolini, D. 2018.", 
                                 "A comparison of bioinformatic approaches for 16S rRNA gene ",
                                 "profiling of food bacterial microbiota. Int. J. Food Microbiol. 265:9-17.",
                                 sep = ""),
                          paste0("Parente, E., De Filippis, F., Ercolini, D., Ricciardi, A., Zotta, T., 2019. ", 
                                 "Advancing integration of data on food microbiome studies: FoodMicrobionet ",
                                 "3.1, a major upgrade of the FoodMicrobionet database. Int. J. Food Microbiol. : 305:108249.",
                                 sep = "")
                          ),
                        DOI = c("10.1016/j.ijfoodmicro.2015.12.001",
                                "10.1016/j.ijfoodmicro.2017.10.028",
                                "10.1016/j.ijfoodmicro.2019.108249")
                        )
saveRDS(references, file = "references.rds")

# app and copyright lines
app_line <- "ShinyFMBN2 (v 2.3) is designed to provide an interface to:"
copyright_line <- "Copyright 2018, 2019, 2020, 2021 Eugenio Parente, Università della Basilicata"


# building and saving the lists -------------------------------------------
# (if all traffic lights are green)

if(any(flag_dupli_sample, flag_dupli_taxa, flag_edges_mismatch_1, 
       flag_edges_mismatch_2, flag_edges_na, flag_sample_mismatch, 
       flag_study_mismatch)){
  cat("WARNING: there are one or more issues with your tables, cannot save the list", "\n")
} else {
  FMBN <- list(
    studies = dplyr::select(studies, studyId:corr_author_mail),
    samples = dplyr::select(samples, studyId:SRA_run),
    edges = edges,
    taxa = taxa,
    references = references,
    version = version_text,
    app_text = app_line,
    copyright_text = copyright_line
  )
  FMBN_plus <- list(
    studies = studies,
    primers = primers,
    samples = samples,
    foodex2exp = foodex2exp,
    edges = edges,
    taxa = taxa,
    references = references,
    version = str_c(version_text, 
                    "this version includes extra fields and tables", 
                    sep = " "),
    app_text = app_line,
    copyright_text = copyright_line
  )
  saveRDS(FMBN, file = "FMBN.rds")
  saveRDS(FMBN_plus, file = "FMBN_plus.rds")
  cat("\nFMBN lists saved", "\n")
}

# Credits and copyright ---------------------------------------------------

# Assume that this is overall under MIT licence

# Copyright 2021 Eugenio Parente
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

