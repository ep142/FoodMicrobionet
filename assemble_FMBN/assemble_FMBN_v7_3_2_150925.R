# assemble_FMBN, a script for assembling multiple files in a working version of FMBN
# LOOK FOR CODE FOR FIXING PROBLEMS IN PREVIOUS VERSIONS OF THE SCRIPT  
# v 7_3_2 15/09/2025
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
# FMBN_F.rds (for fungi) and FMBN_B.rds (for bacteria) which are compatible with ShinyFMBN 2.3 and 3
# FMBN_plus.rds which includes all the data but is not compatible with
# ShinyFMBN 2.3 and 3

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
  cat(red("\nWARNING: One or more .xlsx files are missing!\n"))
}

# version -----------------------------------------------------------------
version_n <- "5.1.2"
last_updated <- "15/09/2025"
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
                     na=c("","NA"))

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
                       na=c("","NA"))


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

# check that all studies are in sample and set a flag
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
  cat(green("\nNo study mismatch detected!\n"))
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
  cat(green("\nNo sample duplication detected!\n"))
}

# check for a rare issue due to wrong copy and paste ops
n_accessions_B <- samples_B %>% 
  group_by(studyId) %>%
  summarise(n_SRR = n_distinct(SRA_run),
            n_biosam = n_distinct(biosample),
            n_SRS = n_distinct(SRA_sample))
n_accessions_B <- left_join(n_accessions_B, select(studies, studyId, samples))
# n_SRR should be equal to samples in study, a mismatch is possible for
# studies which were not downloaded from SRA or in cases in which there
# were errors or incoherences in the data deposited in SRA

saveRDS(samples_B, file = "samples_B.rds")
saveRDS(samples_F, file = "samples_F.rds")


# primers -----------------------------------------------------------------

# use the following command if you want to choose the file manually
# primers_fn <- file.choose()

primers_fn <- "Primers.xlsx"
primers <- read_xlsx(primers_fn,
                     col_types = c(rep("text", 4),
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
  cat(green("\nNo edges with NA values!\n"))
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

# check sums
sum_flags_B <- all(
  near(
    pull(edges_B |> group_by(sampleId) |> summarise(wsum = sum(weight))),
    100,
    tol = 5)
  )
if(!sum_flags_B){
  beep(sound = which_sound)
  cat(red("\nWARNING: For some sampleIds the sum of weights is not close enough to 100!\n"))
} else {
  beep(sound=3)
  cat(green("\nSum of weights in edges_B are within nominal values (ca. 100)!\n"))
}

edges_B <- edges_B |>
  arrange(sampleId, taxonId)

# there is a situation in which 2 edges for the same sample have the same id
# possibly because of fixed SILVA taxonomy
# for brevity will exit after 50 occurrences
flag_dupli_edge_taxId <- F
ndupli <- 0
for(i in seq_along(samples_B$sampleId)){
  mini_edges <- edges_B |>
    dplyr::filter(sampleId == samples_B$sampleId[i])
  if(any(duplicated(mini_edges$taxonId))){
    dupli_taxon <- mini_edges$taxonId[which(duplicated(mini_edges$taxonId))]
    cat("duplicate found", samples_B$sampleId[i], dupli_taxon, "\n")
    ndupli <- ndupli + 1
  }
  if(ndupli >50){
    flag_dupli_edge_taxId <- T
    beep()
    cat(red("\nWARNING: duplicated taxId found for selected sampleId!\n"))
    break
  } 
}

if(flag_dupli_edge_taxId) {
  edges_B_summed <- edges_B |>
    group_by(sampleId, taxonId) |>
    summarise(weight_sum = sum(weight)) |>
    ungroup()
  flag_dupli_edge_taxId <- F
  ndupli <- 0
  for (i in seq_along(samples_B$sampleId)) {
    mini_edges <- edges_B_summed |>
      dplyr::filter(sampleId == samples_B$sampleId[i])
    if (any(duplicated(mini_edges$taxonId))) {
      dupli_taxon <- mini_edges$taxonId[which(duplicated(mini_edges$taxonId))]
      cat("duplicate found",
          samples_B$sampleId[i],
          dupli_taxon,
          "\n")
      if(ndupli >50){
        flag_dupli_edge_taxId <- T
        beep()
        cat(red("\nWARNING: duplicated taxId found for selected sampleId!\n"))
        break
      }
    }
  }
  beep(which_sound)
  cat(green("\nedges_B fixed!\n"))
  
  edges_B_old <- edges_B
  edges_B <- rename(edges_B_summed, weight = weight_sum)
  
  # recheck
  sum_flags_B <- all(
    near(
      pull(edges_B |> group_by(sampleId) |> summarise(wsum = sum(weight))),
      100,
      tol = 5)
  )
  if(!sum_flags_B){
    beep(sound = which_sound)
    cat(red("\nWARNING: For some sampleIds the sum of weights is not close enough to 100!\n"))
  } else {
    beep(sound=3)
    cat(green("\nSum of weights in edges_B are within nominal values (ca. 100)!\n"))
  }
  
}

write_tsv(edges_B, "edges_B_fixed_170925.txt")

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
  cat(green("\nNo mismatch between samples in the samples and edge table for fungi!\n"))
}

# check sums
sum_flags_F <- all(
  near(
    pull(edges_F |> group_by(sampleId) |> summarise(wsum = sum(weight))),
    100,
    tol = 5)
)
if(!sum_flags_F){
  beep(sound = which_sound)
  cat(red("\nWARNING: For some sampleIds the sum of weights is not close enough to 100!\n"))
} else {
  beep(sound=3)
  cat(green("\nSum of weights in edges_F are within nominal values (ca. 100)!\n"))
}

if(!sum_flags_F){
  edges_F_check <- edges_F |>
    summarise(wsum = sum(weight), .by = sampleId)
  edges_F_check <- edges_F_check |>
    dplyr::filter(wsum >101 | wsum <99)
  # duplication possible
  edges_F_to_check <- edges_F |>
    dplyr::filter(sampleId %in% edges_F_check$sampleId) |>
    arrange(sampleId, taxonId)
  # no, other sort of problem, need to check from the original file
  # fix done on 15/10/25 19.15
  edges_minimal_SRP324301 <- read_tsv("edges_minimal_SRP324301_ITS.txt")
  edges_F_temp <- edges_F |>
    dplyr::filter(!(sampleId %in% edges_F_check$sampleId))
  # recheck
  edges_F_check$sampleId %in% edges_F_temp$sampleId
  # add back
  edges_F_temp <- bind_rows(
    edges_F_temp, edges_minimal_SRP324301
  ) |>
    arrange(sampleId, taxonId)
  # recheck
  edges_F_check_2 <- edges_F_temp |>
    summarise(wsum = sum(weight), .by = sampleId)
  edges_F_check_2 <- edges_F_check_2 |>
    dplyr::filter(wsum >101 | wsum <99)
  edges_F <- edges_F_temp
  # recheck
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
    cat(green("\nNo mismatch between samples in the samples and edge table for fungi!\n"))
  }
  
  # check sums
  sum_flags_F <- all(
    near(
      pull(edges_F |> group_by(sampleId) |> summarise(wsum = sum(weight))),
      100,
      tol = 5)
  )
  if(!sum_flags_F){
    beep(sound = which_sound)
    cat(red("\nWARNING: For some sampleIds the sum of weights is not close enough to 100!\n"))
  } else {
    beep(sound=3)
    cat(green("\nSum of weights in edges_F are within nominal values (ca. 100)!\n"))
  }
  
  
}



# there is a situation in which 2 edges for the same sample have the same id
# possibly because of fixed SILVA taxonomy
# for brevity will exit after 50 occurrences
flag_dupli_edge_taxId <- F
ndupli <- 0
for(i in seq_along(samples_F$sampleId)){
  mini_edges <- edges_F |>
    dplyr::filter(sampleId == samples_F$sampleId[i])
  if(any(duplicated(mini_edges$taxonId))){
    dupli_taxon <- mini_edges$taxonId[which(duplicated(mini_edges$taxonId))]
    cat("duplicate found", samples_F$sampleId[i], dupli_taxon, "\n")
    ndupli <- ndupli + 1
  }
  if(ndupli >50){
    flag_dupli_edge_taxId <- T
    beep()
    cat(red("\nWARNING: duplicated taxId found for selected sampleId!\n"))
    break
  } 
}

if(flag_dupli_edge_taxId) {
  edges_F_summed <- edges_F |>
    group_by(sampleId, taxonId) |>
    summarise(weight_sum = sum(weight)) |>
    ungroup()
  flag_dupli_edge_taxId <- F
  ndupli <- 0
  for (i in seq_along(samples_F$sampleId)) {
    mini_edges <- edges_F_summed |>
      dplyr::filter(sampleId == samples_F$sampleId[i])
    if (any(duplicated(mini_edges$taxonId))) {
      dupli_taxon <- mini_edges$taxonId[which(duplicated(mini_edges$taxonId))]
      cat("duplicate found",
          samples_F$sampleId[i],
          dupli_taxon,
          "\n")
      if(ndupli >50){
        flag_dupli_edge_taxId <- T
        beep()
        cat(red("\nWARNING: duplicated taxId found for selected sampleId!\n"))
        break
      }
    }
  }
  beep(which_sound)
  cat(green("\nedges_F fixed!\n"))
  
  edges_F_old <- edges_F
  edges_F <- rename(edges_F_summed, weight = weight_sum)
  
  # recheck
  sum_flags_F <- all(
    near(
      pull(edges_F |> group_by(sampleId) |> summarise(wsum = sum(weight))),
      100,
      tol = 5)
  )
  if(!sum_flags_F){
    beep(sound = which_sound)
    cat(red("\nWARNING: For some sampleIds the sum of weights is not close enough to 100!\n"))
  } else {
    beep(sound=3)
    cat(green("\nSum of weights in edges_F are within nominal values (ca. 100)!\n"))
  }
  
}

write_tsv(edges_F, "edges_F_fixed_150925.txt")
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
  fix_dupli_taxa <- T
} else {
  beep(sound=3)
  cat(green("\nNo duplicated taxa!\n"))
  fix_dupli_taxa <- F
}

flag_dupli_taxonId <- !(sum(duplicated(taxa$taxonId)) == 0)
# which are duplicated?
if(flag_dupli_taxonId){
  dupli_taxa_ids <- taxa$taxonId[duplicated(taxa$taxonId)]
  cat(red("\nWARNING: there are duplicated taxon ids.","\n"))
  cat(red("The duplicated ids are ", dupli_taxa_ids, "\n"))
  cat(red("YOU NEED TO FIX THEM MANUALLY", "\n"))
  fix_dupli_taxonId <- T
} else {
  beep(sound=3)
  cat(green("\nNo duplicated taxon ids!\n"))
  fix_dupli_taxonId <- F
}

if(fix_dupli_taxa){
  taxa_temp <- taxa
  edges_F_temp <- edges_F
  edges_B_temp <- edges_B
  for(i in seq_along(dupli_taxa_labels)){
    dupli_taxon_i <- taxa |> dplyr::filter(label == dupli_taxa_labels[i])
    keep_taxon_i <- dupli_taxon_i$taxonId[1]
    remove_taxon_i <- dupli_taxon_i$taxonId[2]
    taxa_temp <- taxa_temp |> 
      dplyr::filter(taxonId != remove_taxon_i)
    # fix edges
    if(dupli_taxon_i$domain[1]=="Fungi"){
      edges_F_temp <- edges_F_temp  |>
        mutate(taxonId = if_else(taxonId == remove_taxon_i, keep_taxon_i, taxonId))
    } else {
      edges_B_temp <- edges_B_temp  |>
        mutate(taxonId = if_else(taxonId == remove_taxon_i, keep_taxon_i, taxonId))
    }
  }
  taxa <- taxa_temp
  write_tsv(taxa, "taxa_fixed.txt")
  edges_B <- edges_B_temp
  write_tsv(edges_B, "edges_B_fixed.txt")
  edges_F <- edges_F_temp
  write_tsv(edges_F, "edges_F_fixed.txt")
}




# check that taxa match in the edges and taxa table
flag_taxa_B_mismatch <- !(all(unique(edges_B$taxonId) %in%  sort(taxa$taxonId)))
flag_taxa_F_mismatch <- !(all(unique(edges_F$taxonId) %in%  sort(taxa$taxonId)))
refix_B <- F
refix_F <- F
if(any(flag_taxa_B_mismatch, flag_taxa_F_mismatch)){
  if(flag_taxa_B_mismatch){
    refix_B <- T
    cat(red("WARNING: there are taxa in the edge table which are not found in the taxa table", "\n"))
    missing_taxa_B <- which(!(unique(edges_B$taxonId) %in%  sort(taxa$taxonId)))
    cat("The taxonId for the missing bacterial taxa are: ", unique(edges_B$taxonId)[missing_taxa_B])
    if(refix){
      # ad hoc code for fixing edges
      refix_B <- F
    }
    flag_taxa_B_mismatch <- !(all(unique(edges_B$taxonId) %in%  sort(taxa$taxonId)))
    
  }
  if(flag_taxa_F_mismatch){
    refix_F <- T
    cat(red("WARNING: there are taxa in the edge table which are not found in the taxa table", "\n"))
    missing_taxa_F <- which(!(unique(edges_F$taxonId) %in%  sort(taxa$taxonId)))
    cat("\nThe taxonId for the missing fungal taxa are: \n", 
        sort(unique(edges_F$taxonId)[missing_taxa_F]))
    if(refix){
      # ad hoc code
      missing_taxa_lookup_fn <- file.choose()
      # it is an excel file
      missing_taxa_lookup <- read_xlsx(missing_taxa_lookup_fn)
      edges_F_temp <- edges_F
      for(i in seq_along(missing_taxa_lookup$taxon_missing)){
        edges_F_temp <- edges_F_temp |>
          mutate(taxonId = if_else(taxonId == missing_taxa_lookup$taxon_missing[i], 
                                   missing_taxa_lookup$taxon[i], 
                                   taxonId))
        cat("\n replaced taxon ", missing_taxa_lookup$taxon_missing[i], " with ", 
            missing_taxa_lookup$taxon[i])
      }
      !(all(unique(edges_F_temp$taxonId) %in%  sort(taxa$taxonId)))
      edges_F <- edges_F_temp
      write_tsv(edges_F, "edges_F_fixed_taxa.txt")
      refix_F <- F
    }
    flag_taxa_F_mismatch <- !(all(unique(edges_F$taxonId) %in%  sort(taxa$taxonId)))
  } 
} else {
  beep(sound=3)
  cat(green("\nNo missing taxa!\n"))
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
                       paste0("De Filippis, F., Parente, E., Zotta, T., Ercolini, D. 2018. ", 
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
                              "A comprehensive view of food microbiota: introducing FoodMicrobionet v5. ",
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


# app and copyright lines -----------------------------------------------

app_line <- "ShinyFMBN2 (v3) is designed to provide a GUI to:"
copyright_line <- "Copyright 2018, 2019, 2020, 2022, 2023, 2024, 2025 Eugenio Parente, Università della Basilicata; version 5 of FoodMicrobionet was created within the PRIN 2022 project NCYdiversity P20229JMMH, and received funding from the European Union Next-GenerationEU (PIANO NAZIONALE DI RIPRESA E RESILIENZA (PNRR) – MISSIONE 4 COMPONENTE 2, INVESTIMENTO 1.4 – D.D. 1032 17/06/2022). "


# building and saving the list (if all traffic lights are green) --------

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

# Copyright 2021, 2022, 2023, 2024. 2025 Eugenio Parente
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

