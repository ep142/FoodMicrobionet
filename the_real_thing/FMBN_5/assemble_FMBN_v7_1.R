# assemble_FMBN, a script for assembling multiple files in a working version of FMBN

# v 7_1 17/6/2024
# this script takes as an input individual tables for Foodmicrobionet
# Studies
# Primers
# Samples
# FoodEx2exp
# Taxa
# Edges
# Abstracts
# Version history
# The tables must be in .xlsx format (no formulas) and be in the same 
# folder. The exception is the edge table, tab delimited text format
# Creating a project allows to automate loading of the files.
# Version info has to be provided (see line 55)
# Once the tables have been imported in R several diagnostic checks are performed
# to identify potential mismatched in the tables. If the diagnostics fail no
# list is created.
# The tables, together with info on version, are assembled in two types of lists
# FMBN_F.rds (for fungi) and FMBN_B.rds (for bacteria) which are compatible with ShinyFMBN 2.3
# FMBN_plus.rds which includes all the data but is not compatible with
# ShinyFMBN 2.3

# install/load packages

.cran_packages <- c("readxl","tidyverse", "crayon", "beepr")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}

sapply(c(.cran_packages), require, 
       character.only = TRUE)

# set the sound for beepr
which_sound <- 1
# beep(sound = which_sound)

# read filenames --------------------------------------------------------
filenames <- list.files()

# check that files with proper names are available
studies_file_flag <- "Studies.xlsx" %in% filenames
primers_file_flag <- "Primers.xlsx" %in% filenames
samples_b_file_flag <- "Samples_B.xlsx" %in% filenames
samples_f_file_flag <- "Samples_F.xlsx" %in% filenames
FoodEx2_file_flag <- "FoodEx2exp.xlsx" %in% filenames
taxa_file_flag <- "Taxa.xlsx" %in% filenames
edges_file_B_flag <- "Edges_B.txt" %in% filenames
edges_file_F_flag <- "Edges_F.txt" %in% filenames
version_file_flag <- "Version_history.xlsx" %in% filenames
abstract_file_flag <- "Abstracts.xlsx" %in% filenames


file_flags <- c(studies_file_flag, primers_file_flag, samples_b_file_flag,
                samples_f_file_flag, FoodEx2_file_flag, taxa_file_flag, 
                edges_file_B_flag, edges_file_F_flag, version_file_flag,
                abstract_file_flag)

if(all(file_flags)){
  beep(sound = which_sound)
  cat(green("\nAll .xlsx files available!\n"))
} else {
  beep(sound = which_sound)
  cat("\nWARNING: One or more .xlsx files are missing!\n")
}

# version -----------------------------------------------------------------
version_n <- "5.0.1"
last_updated <- "17/06/2024"
version_text <- paste("This is version ", version_n, 
                      " of FoodMicrobionet. Last updated ", last_updated, ".", sep = "")
saveRDS(version_text, file = "version.rds")


# samples -----------------------------------------------------------------

# use the following commands if you want to load the file manually
# samples_fn <- file.choose()
samples_B_fn <- "Samples_B.xlsx"

samples_B <- read_xlsx(samples_B_fn,
                     col_types = c(
                       "text",
                       "numeric",
                       rep("text", 5),
                       rep("numeric", 3),
                       rep("text", 16),
                       "numeric"
                     ),
                     na="NA")

samples_F_fn <- "Samples_F.xlsx"

samples_F <- read_xlsx(samples_F_fn,
                       col_types = c(
                         "text",
                         "numeric",
                         rep("text", 5),
                         rep("numeric", 3),
                         rep("text", 16),
                         "numeric"
                       ),
                       na="NA")


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
              rep("text", 13),
              rep("logical", 2)
            ))

# check if studyId in samples match those in studies
# this must be done for both samples files

studyIds_B <- samples_B |> pull(studyId)
studyIds_F <- samples_F |> pull(studyId)
# union
studyIds_samples <- union(studyIds_B, studyIds_F)

# set a flag
flag_study_mismatch <- !all(sort(studyIds_samples) == sort(studies$studyId))

if(flag_study_mismatch){
  if(length(unique(samples$studyId))>nrow(studies)) {
    beep(sound = which_sound)
    cat(red("\nWARNING: there are more studies in the sample table than in the study table\n"))
    mismatched_studies <- setdiff(unique(samples$studyId), studies$studyId)
    cat(red("\nMismatched studies: ", mismatched_studies,"\n"))
  } else {
    beep(sound = which_sound)
    cat(red("\nWARNING: there are more studies in the study table than in the sample table\n"))
    mismatched_studies <- setdiff(studies$studyId, unique(samples$studyId))
    cat(red("\nMismatched studies: ", mismatched_studies,"\n"))
  }
} else {
  beep(sound=3)
  cat(green("\nNo study mismatch detected\n"))
}

saveRDS(studies, file = "studies.rds")


# check and save samples --------------------------------------------------

# check for duplicates

flag_dupli_sample_B <- !(sum(duplicated(samples_B$label_2))==0)
flag_dupli_sample_F <- !(sum(duplicated(samples_F$label_2))==0)

# find duplicated labels 
if(flag_dupli_sample_B | flag_dupli_sample_F){
  beep(sound = which_sound)
  cat(red("WARNING: there are duplicated samples!!!!\n"))
  if(flag_dupli_sample_B) cat(red("WARNING: the duplicates are in samples for bacteria!!!!\n"))
  if(flag_dupli_sample_F) cat(red("WARNING: the duplicates are in samples for fungi!!!!\n"))
  dupli_sample_labels_B <- samples_B$label_2[duplicated(samples_B$label_2)]
  dupli_sample_labels_F <- samples_F$label_2[duplicated(samples_F$label_2)]
} else {
  beep(sound=3)
  cat(green("\nNo sample duplication detected\n"))
}

# obsolete
# check if sample Id match row numbers
# this is necessary for sample extraction in the Shiny app: mismatches would 
# cause wrong samples to be selected in the aggregate tab

# flag_sample_mismatch <- !(all(samples$sampleId == 1:nrow(samples)))

# if(flag_sample_mismatch) {
#   cat("WARNING: sampleIds DO NOT match row numbers")
#   which(samples$sampleId != 1:nrow(samples))
# }

# check for a rare issue due to wrong copy and paste ops
n_accessions_B <- samples_B %>% 
  group_by(studyId) %>%
  summarise(n_SRR = n_distinct(SRA_run),
            n_biosam = n_distinct(biosample),
            n_SRS = n_distinct(SRA_sample))
n_accessions_B <- left_join(n_accessions_B, select(studies, studyId, samples))
# n_SRR should be equal to samples in study, a mismatch is possible for
# studies which were not downloaded from SRA

saveRDS(samples_B, file = "samples_B.rds")
saveRDS(samples_F, file = "samples_F.rds")


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
# edges_B_fn <- file.choose()
# edges_F_fn <- file.choose()
edges_B_fn <- "Edges_B.txt"
edges_F_fn <- "Edges_F.txt"
# this is the largest file, import might take a while
edges_B <- read_tsv(edges_B_fn)
edges_F <- read_tsv(edges_F_fn)

# check if there are NA values and set a flag
flag_edges_B_na <- anyNA(edges_B)
flag_edges_F_na <- anyNA(edges_F)
if(any(flag_edges_B_na, flag_edges_F_na)) {
  beep(sound = which_sound)
  cat(red("\nWARNING: there are edges with NA values!\n"))
  if(flag_edges_B_na) cat("The edges for bacteria with NA values are at rows ", which(is.na(edges_B)))
  if(flag_edges_F_na) cat("The edges for fungi with NA values are at rows ", which(is.na(edges_F)))
} else {
  beep(sound=3)
  cat(green("\nNo edges with NA values\n"))
}

# edges, bacteria ---------------------------------------------------------

# check if sampleId in edges have the same length as those in samples

# set a flag when number of samples in samples and edges do not match
# for either bacteria or fungi
flag_edges_B_mismatch_1 <- !(length(unique(edges_B$sampleId)) == nrow(samples_B))
# initialize flag for mismatch_2
flag_edges_B_mismatch_2 <- ifelse(flag_edges_B_mismatch_1, T, F)

# check what is wrong
if(flag_edges_B_mismatch_1){
  cat("WARNING: the number of samples in the edges and in the sample table do not match","\n")
  if(length(unique(edges_B$sampleId))> nrow(samples_B)){
    cat("the excess edges are in the edge table at rows: ", setdiff(unique(edges_B$sampleId), samples_B$sampleId))
  } else {
    cat("the missing edges are in the sample table, samples: ", setdiff(samples_B$sampleId, unique(edges_B$sampleId)))
  }
} else {
  # check that sampleId are identical in both tables
  all(sort(unique(edges_B$sampleId)) == sort(samples_B$sampleId))
  flag_edges_mismatch_2 <- !all(sort(unique(edges_B$sampleId)) == sort(samples_B$sampleId))
  if(flag_edges_mismatch_2) cat("the mismatching edges are ", setdiff(unique(edges_B$sampleId), samples_B$sampleId))
}

if(!flag_edges_B_mismatch_1){
  beep(sound=3)
  cat(green("\nNo mismatch between samples in the samples and edge table for bacteria\n"))
}

run_me <- F # this was fixed

if(run_me) {
  # edges mismatch, the edges from studies 236 (ERP146095) and 250 (SRP346843) are missing
  # fixing
  # opening the edges for ST136
  st136_edges_fn <- file.choose()
  st136_edge <- read_tsv(st136_edges_fn)
  edges_B_v2 <-
    bind_rows(edges_B, select(st136_edge, sampleId, taxonId, weight)) |>
    arrange(sampleId, taxonId)
  st150_edges_fn <- file.choose()
  st150_edge <- read_tsv(st150_edges_fn)
  edges_B_v2 <-
    bind_rows(edges_B_v2, select(st150_edge, sampleId, taxonId, weight)) |>
    arrange(sampleId, taxonId)
  # recheck
  all(samples_B$sampleId == unique(edges_B_v2$sampleId))
  length(samples_B$sampleId)
  length(unique(edges_B_v2$sampleId))
  # there are edges still missing
  # for either bacteria or fungi
  flag_edges_B_mismatch_1 <-
    !(length(unique(edges_B_v2$sampleId)) == nrow(samples_B))
  # initialize flag for mismatch_2
  flag_edges_B_mismatch_2 <- ifelse(flag_edges_B_mismatch_1, T, F)
  
  # check what is wrong
  if (flag_edges_B_mismatch_1) {
    cat("WARNING: the number of samples in the edges and in the sample table do not match",
        "\n")
    if (length(unique(edges_B_v2$sampleId)) > nrow(samples_B)) {
      cat(
        "the excess edges are in the edge table at rows: ",
        setdiff(unique(edges_B_v2$sampleId), samples_B$sampleId)
      )
    } else {
      cat(
        "the missing edges are in the sample table, samples: ",
        setdiff(samples_B$sampleId, unique(edges_B_v2$sampleId))
      )
    }
  } else {
    # check that sampleId are identical in both tables
    all(sort(unique(edges_B_v2$sampleId)) == sort(samples_B$sampleId))
    flag_edges_B_mismatch_2 <-
      !all(sort(unique(edges_B_v2$sampleId)) == sort(samples_B$sampleId))
    if (flag_edges_B_mismatch_2)
      cat("the mismatching edges are ",
          setdiff(unique(edges_B_v2$sampleId), samples_B$sampleId))
  }
  
  edges_B <- edges_B_v2
  # triple check
  all(samples_B$sampleId == unique(edges_B$sampleId))
  if(!all(samples_B$sampleId == unique(edges_B$sampleId))){
    beep(sound=which_sound)
    cat(red("\nWARNING: there is still a mismatch between samples in the sample table and samples in the edge table"))
  } else {
    beep(sound=3)
    cat(green("\nNo mismatch between samples in the samples and edge table for bacteria\n"))
  }
  write_tsv(edges_B, file = "FMBN_5_edges_B.txt")
  write_tsv(edges_B, file = "Edges_B.txt")
}



saveRDS(edges_B, file = "edges_B.rds")


# edges, fungi ------------------------------------------------------------


# check if sampleId in edges have the same length as those in samples

# set a flag when number of samples in samples and edges do not match
# for either bacteria or fungi
flag_edges_F_mismatch_1 <- !(length(unique(edges_F$sampleId)) == nrow(samples_F))
# initialize flag for mismatch_2
flag_edges_F_mismatch_2 <- ifelse(flag_edges_F_mismatch_1, T, F)

# check what is wrong
if(flag_edges_F_mismatch_1){
  cat(red("\nWARNING: the number of samples in the edges and in the sample table do not match","\n"))
  if(length(unique(edges_F$sampleId))> nrow(samples_F)){
    cat("the excess edges are in the edge table at rows: ", 
        setdiff(unique(edges_F$sampleId), samples_F$sampleId))
  } else {
    cat("the missing edges are in the sample table, samples: ",
        setdiff(samples_F$sampleId, unique(edges_F$sampleId)))
  }
} else {
  # check that sampleId are identical in both tables
  all(sort(unique(edges_F$sampleId)) == sort(samples_F$sampleId))
  flag_edges_F_mismatch_2 <- !all(sort(unique(edges_F$sampleId)) == sort(samples_F$sampleId))
  if(flag_edges_F_mismatch_2) cat(red("the mismatching edges are ", setdiff(unique(edges_F$sampleId, samples_F$sampleId))))
}

if(!flag_edges_F_mismatch_1){
  beep(sound=3)
  cat(green("\nNo mismatch between samples in the samples and edge table for bacteria\n"))
}


# the missing edges are in study ST234, SRP400743
run_me_F <- F # this was fixed

if(run_me_F) {
  # edges mismatch, the edges from studies 234 (SRP400743) are missing
  # fixing
  # opening the edges for ST136
  st234_edges_fn <- file.choose()
  st234_edge <- read_tsv(st234_edges_fn)
  edges_F_v2 <-
    bind_rows(edges_F, select(st234_edge, sampleId, taxonId, weight)) |>
    arrange(sampleId, taxonId)
  # recheck
  all(samples_F$sampleId == unique(edges_F_v2$sampleId))
  length(samples_F$sampleId)
  length(unique(edges_F_v2$sampleId))

  flag_edges_F_mismatch_1 <-
    !(length(unique(edges_F_v2$sampleId)) == nrow(samples_F))
  # initialize flag for mismatch_2
  flag_edges_F_mismatch_2 <- ifelse(flag_edges_F_mismatch_1, T, F)
  
  # check what is wrong
  if (flag_edges_F_mismatch_1) {
    cat("WARNING: the number of samples in the edges and in the sample table do not match",
        "\n")
    if (length(unique(edges_F_v2$sampleId)) > nrow(samples_F)) {
      cat(
        "the excess edges are in the edge table at rows: ",
        setdiff(unique(edges_F_v2$sampleId), samples_F$sampleId)
      )
    } else {
      cat(
        "the missing edges are in the sample table, samples: ",
        setdiff(samples_F$sampleId, unique(edges_F_v2$sampleId))
      )
    }
  } else {
    # check that sampleId are identical in both tables
    all(sort(unique(edges_F_v2$sampleId)) == sort(samples_F$sampleId))
    flag_edges_mismatch_2 <-
      !all(sort(unique(edges_F_v2$sampleId)) == sort(samples_F$sampleId))
    if (flag_edges_mismatch_2)
      cat("the mismatching edges are ",
          setdiff(unique(edges_F_v2$sampleId), samples_F$sampleId))
  }
  
  edges_F <- edges_F_v2
  # triple check
  all(samples_F$sampleId == unique(edges_F$sampleId))
  if(!all(samples_F$sampleId == unique(edges_F$sampleId))){
    beep(sound=which_sound)
    cat(red("\nWARNING: there is still a mismatch between samples in the sample table and samples in the edge table"))
  } else {
    beep(sound=3)
    cat(green("\nNo mismatch between samples in the samples and edge table for bacteria\n"))
  }
  write_tsv(edges_F, file = "FMBN_5_edges_F.txt")
  write_tsv(edges_F, file = "Edges_F.txt")
}


saveRDS(edges_F, file = "edges_F.rds")

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
  cat(red("\nWARNING: there are duplicated taxa labels.","\n"))
  cat(red("The duplicated taxa are ", dupli_taxa_labels, "\n"))
} else {
  beep(sound=3)
  cat(green("\nNo duplicated taxa\n"))
}


# check that taxa match in the edges and taxa table
flag_taxa_B_mismatch <- !(all(unique(edges_B$taxonId) %in%  sort(taxa$taxonId)))
flag_taxa_F_mismatch <- !(all(unique(edges_F$taxonId) %in%  sort(taxa$taxonId)))
if(any(flag_taxa_B_mismatch,flag_taxa_F_mismatch)){
  cat("WARNING: there are taxa in the edge table which are not found in the taxa table", "\n")
  missing_taxa_B <- which(!(unique(edges_B$taxonId) %in%  sort(taxa$taxonId)))
  missing_taxa_F <- which(!(unique(edges_F$taxonId) %in%  sort(taxa$taxonId)))
  cat("The taxonId for the missing bacterial taxa are: ", unique(edges_B$taxonId)[missing_taxa_B])
  missing_taxa_F <- which(!(unique(edges_F$taxonId) %in%  sort(taxa$taxonId)))
  cat("The taxonId for the missing fungal taxa are: ", unique(edges_F$taxonId)[missing_taxa_F])
} else {
  beep(sound=3)
  cat(green("\nNo missing taxa\n"))
}


saveRDS(taxa, file = "taxa.rds")


# version hystory ---------------------------------------------------------

version_history_fn <- "Version_history.xlsx"
version_history <- read_xlsx(version_history_fn)
saveRDS(version, file = "version_history.rds")


# abstract --------------------------------------------------------------

abstracts_fn <- "Abstracts.xlsx"
abstracts <- read_xlsx(abstracts_fn)
saveRDS(abstracts, file = "abstracts.rds")


# references --------------------------------------------------------------

references <- tibble(studyId = rep("Foodmicrobionet",5),
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
                              "3.1, a major upgrade of the FoodMicrobionet database. Int. J. Food Microbiol. 305:108249.",
                              sep = ""),
                       paste0("Parente, E., Zotta, T., Ricciardi, A., 2022. ", 
                              "FoodMicrobionet v4: A large, integrated, open and transparent database for food bacterial communities. ",
                              "Int J Food Microbiol 372, 109696.",
                              sep = ""),
                       paste0("Parente, E., Ricciardi, A., 2024. ", 
                              "A comprehensive view of food microbiota: introducing FoodMicrobionet v5",
                              "Foods, 13, 1689.",
                              sep = "")
                       ),
                     DOI = c("10.1016/j.ijfoodmicro.2015.12.001",
                             "10.1016/j.ijfoodmicro.2017.10.028",
                             "10.1016/j.ijfoodmicro.2019.108249",
                             "10.1016/j.ijfoodmicro.2022.109696",
                             "10.3390/foods1311168")
)
saveRDS(references, file = "references.rds")


# app and copyright lines
app_line <- "ShinyFMBN2 (v 2.5) is designed to provide an interface to:"
copyright_line <- "Copyright 2018, 2019, 2020, 2022, 2023, 2024 Eugenio Parente, Università della Basilicata; version 5 of FoodMicrobionet was created within the PRIN 2022 project NCYdiversity P20229JMMH, and received funding from the European Union Next-GenerationEU (PIANO NAZIONALE DI RIPRESA E RESILIENZA (PNRR) – MISSIONE 4 COMPONENTE 2, INVESTIMENTO 1.4 – D.D. 1032 17/06/2022). "

# building and saving the list (if all traffic lights are green)

if(any(flag_dupli_sample_B, flag_dupli_sample_F, flag_dupli_taxa, 
       flag_edges_B_mismatch_1, flag_edges_B_mismatch_2, 
       flag_edges_F_mismatch_1, flag_edges_F_mismatch_2, 
       flag_edges_B_na, flag_edges_F_na,flag_study_mismatch,
       flag_taxa_B_mismatch, flag_taxa_F_mismatch)){
  cat(red("WARNING: there are one or more issues with your tables, cannot save the list", "\n"))
} else {
  beep(sound=3)
  cat(green("\nYour list is good to go!\n"))
  studies_F <- studies |>
    dplyr::filter(str_detect(tax_database, "UNITE"))
  FMBN_F <- list(
    studies = dplyr::select(studies_F, studyId:corr_author_mail),
    samples = dplyr::select(samples_F, studyId:SRA_run),
    edges = edges_F,
    taxa = taxa,
    references = references,
    version = version_text,
    app_text = app_line,
    copyright_text = copyright_line
  )
  studies_B <- studies |>
    dplyr::filter(str_detect(target, "16S"))
  FMBN_B <- list(
    studies = dplyr::select(studies_B, studyId:corr_author_mail),
    samples = dplyr::select(samples_B, studyId:SRA_run),
    edges = edges_B,
    taxa = taxa,
    references = references,
    version = version_text,
    app_text = app_line,
    copyright_text = copyright_line
  )
  FMBN_plus <- list(
    studies = studies,
    primers = primers,
    samples_B = samples_B,
    samples_F = samples_F,
    foodex2exp = foodex2exp,
    edges_B = edges_B,
    edges_F = edges_F,
    taxa = taxa,
    abstracts = abstracts,
    version_history = version_history,
    references = references,
    version = str_c(version_text, 
                    "this version includes extra fields and tables", 
                    sep = " "),
    app_text = app_line,
    copyright_text = copyright_line
  )
  saveRDS(FMBN_B, file = "FMBN_B.rds")
  cat("FMBN_B list saved", "\n")
  saveRDS(FMBN_F, file = "FMBN_F.rds")
  cat("FMBN_F list saved", "\n")
  saveRDS(FMBN_plus, file = "FMBN_plus.rds")
  cat("FMBN_plus list saved", "\n")
}


# Credits and copyright ---------------------------------------------------

# Assume that this is overall under MIT licence

# Copyright 2021, 2022, 2023, 2024 Eugenio Parente
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

# version 5 of FoodMicrobionet was created within the PRIN 2022 project 
# NCYdiversity P20229JMMH, and received funding from the European Union 
# Next-GenerationEU (PIANO NAZIONALE DI RIPRESA E RESILIENZA (PNRR) – 
# MISSIONE 4 COMPONENTE 2, INVESTIMENTO 1.4 – D.D. 1032 17/06/2022).

