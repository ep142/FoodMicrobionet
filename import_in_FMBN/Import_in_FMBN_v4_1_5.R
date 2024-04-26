# Import_in_FMBN ----------------------------------------------------------

# a script for assisting the inclusion of new studies in FMBN
# version 4.1.5 23/3/2024

# WARNING the script will not work with R < 4.1


# loading packages --------------------------------------------------------

cran_packages <- c("readxl", "tidyverse", "crayon", "magrittr", "beepr")
sapply(cran_packages, require, character.only = T)

which_sound <- 1

# warning and checking files ----------------------------------------------
beep(which_sound)
cat(red("\nWARNING\nBEFORE YOU BEGIN: copy the mindata file to the appropriate location\n"))
beep(which_sound)
cat(red("\nWARNING\nBEFORE YOU BEGIN: if you are working on bacteria copy the support folder in the project folder\n"))

# all the files need to be in the same (project) folder
# get file list
file_list <- list.files()
# load the edges
edge_file <- file_list[!is.na(str_match(file_list, "edge_table"))]
edges <- read_tsv(edge_file)
# number of samples in edge table (for cross checking)
n_samples_edge <- length(unique(edges$Run_s))
# number of taxa in edge table (for cross checking)
n_taxa_edge <- length(unique(edges$s_label))

beep(which_sound)
cat(red("\nWhat data are you importing? (B for bacteria, F for fungi)\n"))
what_data <- "B"

# Studies -----------------------------------------------------------------

# import the study file

bioproject <- "PRJNA781236"

study_file <- file_list[!is.na(str_match(file_list, "_study.txt"))]

# open, reorder the columns and save

study <- read_tsv(study_file)
if("tax_database" %in% colnames(study)) {
  study <- select(study, - tax_database)}

beep(which_sound)
cat(red("\nWARNING\ncheck if study", study$Seq_accn, "is in the study table and set in_FMBN accordingly\n"))

in_FMBN <- F

if(in_FMBN){
  # need to check if the same accession was used for bioproject
  beep(which_sound)
  cat(red("\nWARNING\nManually edit your main study file: region, primers, accessions..."))
}

# ad hoc
# study$target <- "ITS region"


beep(which_sound)
cat(red("\nWARNING:\nthis condition must be true\n"))
# cross-check number of samples, expression must be TRUE
study$samples == n_samples_edge

if(!("overlapping" %in% colnames(study))) study$overlapping <- NA
if(!("paired_end" %in% colnames(study))) study$paired_end <- NA

# make pipeline fields

# tax_database must be set manually
# if the study includes only bacteria use SILVA v138.1 only
# if the study includes only fungi use UNITEgr250723 only
# otherwise use both "SILVA v138.1 + UNITEgr250723"

pip_fields <- tibble(bioinf_software = "R dada2",
                     OTU_picking = "ASV with DADA2",
                     assign_tax_method = "assignTaxonomy()",
                     tax_database = "SILVA v138.1")
empty_fields <- tibble(study = NA, 
                       studyId = NA, 
                       FMBN_version = NA,
                       Seq_accn_link = NA,
                       bioproject = bioproject,
                       food_group = NA,
                       short_descr = NA,
                       DOI_link = NA,
                       ref_short = NA,
                       year = NA,
                       ref_complete = NA,
                       corr_author_surname = NA, 
                       corr_author_mail = NA,
                       study_type = NA)


study <- bind_cols(study, pip_fields, empty_fields) |>
  select(study,	studyId, FMBN_version, target, region, platform,	
         read_length_bp, seq_center,	bioinf_software,	OTU_picking,	
         assign_tax_method,	tax_database,	Seq_accn,	Seq_accn_link,	bioproject,	
         samples,	food_group, short_descr, DOI, DOI_link, ref_short, year, 
         ref_complete, corr_author_surname, corr_author_mail, geoloc, study_type,
         primer_f, primer_r, overlapping, paired_end)

# export the file
if(!in_FMBN){
  write_tsv(study, paste(study$Seq_accn, "_study_FMBN.txt", sep = ""))
  }

beep(which_sound)
cat(red("\nWARNING\nCheck primers and geolocation and manually edit missing fields; add the abstract to the Abstracts table, if needed\n"))   

# samples -----------------------------------------------------------------
beep(which_sound)
cat(red("\nWARNING\nCheck if the same study is already available in FMBN and if some or all samples are there; 
if they are, match by biosample, if at all possible"))
        

# from Excel files FMBNxxx_studies, FMBNxxx_samples
studyId <- "ST250"
max_samples <- 14929
first_sample <- max_samples +1
# max_samples must be the highest between the maximum number of samples if the
# FMBN_samples_B and FMBN_sample_F

samples_file <- file_list[!is.na(str_match(file_list, "_samples.txt"))]
samples <- read_tsv(samples_file)

n_samples <- nrow(samples)
rep_NA <- rep(NA, n_samples)
rep_NA_integer <- rep(NA_integer_, n_samples)
# double check samples number
study$samples == n_samples # must be T

if(!in_FMBN){
  missing_col_samples <- tibble(studyId = rep(studyId, n_samples),
                                sampleId = seq(from = first_sample,
                                               to = (first_sample+n_samples-1)),
                                label_1 = rep_NA,
                                label_3 = rep_NA,
                                llabel = rep_NA,
                                s_type = rep_NA,
                                n_reads = rep_NA,
                                foodId = rep_NA,
                                L1 = rep_NA,
                                L4 = rep_NA,
                                L6 = rep_NA,
                                nature = rep_NA,
                                process = rep_NA,
                                spoilage = rep_NA,
                                matchId = rep_NA_integer
  )
  
  # this needs to be adjusted manually, to double check inconsistencies
  samples <- bind_cols(samples, missing_col_samples) |>
    mutate(n_reads = n_reads2) |>
    select(studyId,	sampleId,	label_1,	label_2 = label2,	label_3,	llabel,	s_type,	
           n_reads, n_reads2,	n_issues, foodId,	description,	L1,	L4,	L6,	nature,	
           process, spoilage,	target1,	target2,	biosample,	
           SRA_Sample,	SRA_run, geo_loc_country, geo_loc_continent, lat_lon, matchId)
} else {
  if(what_data == "B") {
    beep(which_sound)
    cat(red("\nManually find the position of the main FMBN_samples_F file\n"))
  } else {
    beep(which_sound)
    cat(red("\nManually find the position of the main FMBN_samples_B file\n"))  
  }
  FMBN_samples_fn <- file.choose()
  FMBN_samples <- read_xlsx(FMBN_samples_fn, sheet = "samples")
  # match the samples
  # the samples which are already there
  biosamples_in_dataset <- pull(samples, biosample)
  samples_found <- any(biosamples_in_dataset %in% FMBN_samples$biosample)
  if(!samples_found) {
    beep(which_sound)
    cat(red("Could not find matching biosamples, creating new samples"))
    missing_col_samples <- tibble(studyId = rep(studyId, n_samples),
                                  sampleId = seq(from = first_sample,
                                                 to = (first_sample+n_samples-1)),
                                  label_1 = rep_NA,
                                  label_3 = rep_NA,
                                  llabel = rep_NA,
                                  s_type = rep_NA,
                                  n_reads = rep_NA,
                                  foodId = rep_NA,
                                  L1 = rep_NA,
                                  L4 = rep_NA,
                                  L6 = rep_NA,
                                  nature = rep_NA,
                                  process = rep_NA,
                                  spoilage = rep_NA,
                                  matchId = rep_NA_integer
    )
    
    # this needs to be adjusted manually, to double check inconsistencies
    samples <- bind_cols(samples, missing_col_samples) |>
      mutate(n_reads = n_reads2) |>
      select(studyId,	sampleId,	label_1,	label_2 = label2,	label_3,	llabel,	s_type,	
             n_reads, n_reads2,	n_issues, foodId,	description,	L1,	L4,	L6,	nature,	
             process, spoilage,	target1,	target2,	biosample,	
             SRA_Sample,	SRA_run, geo_loc_country, geo_loc_continent, lat_lon,
             matchId)  
  } else {
    # filter the FMBN file
    FMBN_samples_in_dataset <- FMBN_samples |>
      dplyr::filter(biosample %in% biosamples_in_dataset)
    # remove unnecessary info
    FMBN_samples_in_dataset <- FMBN_samples_in_dataset |>
      dplyr::select(biosample, studyId, sampleId, label_1, label_3, llabel, s_type, foodId, 
                    description, L1, L4, L6, nature, process, spoilage)
    
    
    warn_dupli_biosample <- if_else(
      length(unique(FMBN_samples_in_dataset$biosample)) < nrow(FMBN_samples_in_dataset) |
        length(unique(samples$biosample)) < nrow(samples),
      T, F
    )
    if(warn_dupli_biosample){
      beep(which_sound)
      cat(red("\n***WARNING:\n there are duplicated biosample accessions in one of your sample files.\n
              You need to write the code for sample matching manually, NO SAMPLE TABLE WILL BE SAVED OTHERWISE"))
      
      # WRITE THE CODE FOR SAMPLE MATCHING HERE
    } else {
      # join the existing sample info based on biosample
      samples_plus <- left_join(rename(samples, description_2 = description), 
                                FMBN_samples_in_dataset)
      samples_plus <- samples_plus |>
        mutate(n_reads = n_reads2) |>
        select(studyId,	sampleId,	label_1,	label_2 = label2,	label_3,	llabel,	s_type,	
               n_reads, n_reads2,	n_issues, foodId,	description,	L1,	L4,	L6,	nature,	
               process, spoilage,	target1,	target2,	biosample,	
               SRA_Sample,	SRA_run, geo_loc_country, geo_loc_continent, lat_lon)
      
      # any sample not in FMBN?
      if(nrow(FMBN_samples_in_dataset) < nrow(samples_plus)){
        beep(which_sound)
        cat(red("\nWARNING: you have more samples than there are in FMBN, will add sample Ids"))
        
        excess_samples <- nrow(samples_plus)-nrow(FMBN_samples_in_dataset)
        
        beep(which_sound)
        cat(red("\nWARNING: need to add ", excess_samples, " to the study entry in FMBN"))
        cat(red("\nthe number of total samples in the study is actually",
                length(union(FMBN_samples_in_dataset$biosample, samples$biosample)),"\n"))

        samples_to_Id <- samples_plus |> dplyr::filter(is.na(sampleId)) 
        samples_to_Id$sampleId <- seq(from = first_sample, 
                                      to = (first_sample + (nrow(samples_to_Id)-1)))
        samples_to_Id$studyId <- studyId
        descriptions <- samples |>
          dplyr::filter(SRA_run %in% samples_to_Id$SRA_run) |>
          arrange(SRA_run)
        samples_to_Id$description <- descriptions$description
        samples_plus <- dplyr::filter(samples_plus, !is.na(sampleId))
        samples_plus <- bind_rows(samples_plus, samples_to_Id) |>
          arrange(SRA_run)
      }
      # add matchId 
      samples_plus <- left_join(samples_plus, 
                                select(FMBN_samples_in_dataset, biosample, matchId = sampleId))
      
      # build the id matches to put in the sample file of B (if F) or F (if B)
      matchedIds <- samples_plus |>
        dplyr::filter(!is.na(matchId)) |>
        mutate(sampleId_new = sampleId) |>
        select(-sampleId) |>
        dplyr::select(sampleId = matchId, match_to_old =sampleId_new) |>
        arrange(sampleId)
      write_tsv(matchedIds, paste(study$Seq_accn, "_match_samples_FMBN.txt", sep = ""))
      samples_plus <- samples_plus |>
        select(studyId,	sampleId,	label_1,	label_2,	label_3,	llabel,	s_type,	
               n_reads, n_reads2,	n_issues, foodId,	description,	L1,	L4,	L6,	nature,	
               process, spoilage,	target1,	target2,	biosample,	
               SRA_Sample,	SRA_run, geo_loc_country, geo_loc_continent, lat_lon, matchId)
      samples <- samples_plus
      }
  }
}

# export the file
write_tsv(samples, paste(study$Seq_accn, "_samples_FMBN.txt", sep = ""))

# taxa --------------------------------------------------------------------
beep(which_sound)
cat(red("\nWARNING:\ncreate an Excel file with taxa from the taxa table.\n"))

taxa_in_FMBN <- read_excel("taxa_in_fmbn.xlsx")

max_taxaId <- max(taxa_in_FMBN$taxonId)

beep(which_sound)
cat(red("\nWARNING:\ncheck that there are no duplicates: must return T and 0.\n"))
length(unique(taxa_in_FMBN$label))==length(taxa_in_FMBN$label)
anyDuplicated(taxa_in_FMBN$label)

# remove duplicates
if(length(taxa_in_FMBN$label[anyDuplicated(taxa_in_FMBN$label)])>0){
  print("The following are duplicated and will be removed:")
  taxa_in_FMBN$label[anyDuplicated(taxa_in_FMBN$label)]
  taxa_in_FMBN <- taxa_in_FMBN[-anyDuplicated(taxa_in_FMBN$label),]
}


# open taxa

taxa_file <- file_list[!is.na(str_match(file_list, "taxaFMBN"))]
taxa <- read_tsv(taxa_file)

# do cross checks, all must be TRUE
beep(which_sound)
cat(red("\nWARNING:\nthis condition must be true\n"))
nrow(taxa) == n_taxa_edge

# if it fails at first check, need to recheck the original edge file
# still fails, must see which is missing in the edges
length(unique(edges$s_label))
which(duplicated(taxa$label))
# duplication is generally due to "Incertae Sedis" in genus and label
if(length(which(taxa$Genus == "Incertae Sedis"))>0){
  pos_to_change <- which(taxa$Genus == "Incertae Sedis")
  taxa[pos_to_change, ] <- taxa |> dplyr::filter(Genus == "Incertae Sedis") |>
    mutate(label = str_c(Family, Genus, sep = " ")) |>
    mutate(Genus = str_c(Family, Genus, sep = " ")) |>
    mutate(taxonomy = str_c("k__", Kingdom, "; p__", Phylum, "; c__", Class,
                            "; o__; f__", Family, "; g__", Genus,
                            "; s__")) |>
    mutate(taxonomy_L6 = str_c("k__", Kingdom, "; p__", Phylum, "; c__", Class,
                            "; o__", Order, "; f__", Family, "; g__", Genus,
                            "; s__")) |>
    mutate(id = str_c("Root;k__", Kingdom, "__", Phylum, ";c__", Class,
                            ";f__", Family, ";g__", Genus,
                            ";s__")) |>
    mutate(id_L6 = str_c("Root;k__", Kingdom, "__", Phylum, ";c__", Class,
                      ";o__", Order, ";f__", Family, ";g__", Genus,
                      ";s__"))
}
# you may also try (after visually inspecting)
# taxa <- taxa[-which(duplicated(taxa$label)),]
# which(duplicated(taxa$label))

all(sort(unique(edges$s_label)) == sort(taxa$label))
 
# check this if it fails
# fails, even if same length
# check_taxa <- tibble(edges = sort(unique(edges$s_label)),
#                      taxa = sort(taxa$label))
# again a problem with Actinobacteria
# taxa[taxa$label=="Actinobacteria",]$label <- "Actinobacteria (class)"
# this is temporary, I need to fix it before
if(length(which(edges$s_label == "Incertae Sedis"))>0){
  pos_to_change_e <- which(edges$s_label == "Incertae Sedis")
  edges$s_label[pos_to_change_e] <- taxa$Genus[pos_to_change[1]]
}

# recheck
check_taxa <- tibble(edges = sort(unique(edges$s_label)),
                     taxa = sort(taxa$label))

check_taxa <- check_taxa |> 
  mutate(issue= (edges != taxa))

beep(which_sound)
cat(red("\nWARNING:\nthis condition must be true\n"))
all(sort(unique(edges$s_label)) == sort(taxa$label))

# change selected taxa and edges ------------------------------------------

# fix a few taxa (need to be updated every time I spot a new one)
# preserve original versions of taxa and edges
taxa_original <- taxa
edges_original <- edges
# check that all weight sums are 100, must be true
beep(which_sound)
cat(red("\nWARNING:\nthis condition must be true\n"))
all(near(pull(edges_original |> group_by(Run_s) |> summarise(wsum = sum(weight))), 100))

if("mitochondria" %in% taxa$label){
  to_change <- which(taxa$label == "mitochondria")
  taxa$label[to_change,] <- "Mitochondria"
  to_change <- which(edges$s_label == "mitochondria")
  edges$s_label[to_change] <- "Mitochondria"
  print("mitochondria changed")
}

if(what_data == "B"){
  if("Clostridium_sensu_stricto 1" %in% taxa$Genus){
    to_change <- which(taxa$Genus == "Clostridium_sensu_stricto 1")
    taxa$Genus[to_change] <- "Clostridium_sensu_stricto_1"
    taxa[to_change] |>
      mutate(label = if_else(label == "Clostridium_sensu_stricto 1", 
                             "Clostridium_sensu_stricto_1",
                             label))
    to_change <- which(edges$s_label == "Clostridium_sensu_stricto 1")
    edges$s_label[to_change] <- "Clostridium_sensu_stricto_1"
    print("Clostridium_sensu_stricto 1 changed")
  }
  
  if("Clostridium_sensu_stricto 12" %in% taxa$Genus){
    to_change <- which(taxa$Genus == "Clostridium_sensu_stricto 12")
    taxa$Genus[to_change] <- "Clostridium_sensu_stricto_12"
    taxa[to_change] |>
      mutate(label = if_else(label == "Clostridium_sensu_stricto 12", 
                             "Clostridium_sensu_stricto_12",
                             label))
    to_change <- which(edges$s_label == "Clostridium_sensu_stricto 12")
    edges$s_label[to_change] <- "Clostridium_sensu_stricto_12"
    print("Clostridium_sensu_stricto 12 changed")
  }
  
  if("Escherichia-Shigella" %in% taxa$Genus){
    to_change <- which(taxa$Genus == "Escherichia-Shigella")
    taxa$Genus[to_change] <- "Escherichia/Shigella"
    taxa[to_change,] %<>%
      mutate(label = if_else(label == "Escherichia-Shigella", 
                             "Escherichia/Shigella",
                             label))
    to_change <- which(edges$s_label == "Escherichia-Shigella")
    edges$s_label[to_change] <- "Escherichia/Shigella"
    print("Escherichia-Shigella changed")
  }
  
  if("Incertae_Sedis" %in% taxa$label){
    to_change <- which(taxa$label == "Incertae_Sedis" & taxa$Family == "Ruminococcaceae")
    taxa$Genus[to_change] <- NA_character_
    taxa$label[to_change] <- "Ruminococcaceae Incertae Sedis"
    to_change <- which(edges$s_label == "Incertae_Sedis")
    edges$s_label[to_change] <- "Ruminococcaceae Incertae Sedis"
    print("Incertae sedis changed")
  }
  
  # handle more changes using a lookup table
  
  lookup_table <- read_tsv(file.path("support","species_lookup.txt"))
  colnames(lookup_table)[1]<-"label"
  to_change <- which(lookup_table$label %in% taxa$label)
  
  taxa <- left_join(taxa, lookup_table)
  taxa <- taxa |>
    mutate(label = ifelse(is.na(new_id), label, new_id),
           Species = ifelse(is.na(new_species), Species, new_species),
           Genus = ifelse(is.na(new_genus), Genus, new_genus),
           Family = ifelse(is.na(new_family), Family, new_family)
    ) |>
    select(-(new_id:new_family))
  
  edges <- left_join(edges, select(lookup_table, label, new_id), 
                     by=c("s_label" = "label"))
  edges <- edges |>
    mutate(s_label = ifelse(is.na(new_id), s_label, new_id)) |>
    select(-new_id)
  
  # report changes
  if(length(to_change)>0){
    cat("The following taxa were changed to make them coherent with FMBN","\n")
    lookup_table$label[to_change]
  } else {
    cat("No changes were made to taxa using the lookup table")
  }
}

# taxa not in FMBN --------------------------------------------------------

taxa_not_in_FMBN <- anti_join(taxa, select(taxa_in_FMBN, label))

# taxa with changed lineages in SILVA 138 or UNITE
taxa_ch_lineage <- inner_join(select(taxa, label:Species, id_L6), select(taxa_in_FMBN, label, id_L6_old = id_L6)) |>
  dplyr::filter(id_L6 != id_L6_old) |>
  dplyr::filter(!(label %in% taxa_not_in_FMBN$label))

# double check
in_both <- intersect(taxa_in_FMBN$label, taxa$label)
in_taxa <- setdiff(taxa$label, taxa_in_FMBN$label)

# must be true 
beep(which_sound)
cat(red("\nWARNING:\nthis condition must be true\n"))
length(in_both) + length(in_taxa) == nrow(taxa)

# recheck
dupli_edges <- F
if(length(unique(taxa$label))<length(taxa$label)){
  dupli <- duplicated(taxa$label)
  cat("The following are duplicated in taxa", taxa$label[taxa$label %in% taxa$label[dupli]],
      "and will be removed")
  # there might be a need for checking edges for duplicates
  dupli_edges<-T
  # now remove duplicates from taxa
  taxa <- taxa[-anyDuplicated(taxa$label),]
}

if (nrow(taxa_not_in_FMBN)>0){
  taxa_to_FMBN <- taxa_not_in_FMBN |>
    mutate(taxonId = seq(from = max_taxaId+1, 
                         to = max_taxaId+nrow(taxa_not_in_FMBN))) |>
    select(taxonId,	label,	domain = Kingdom,	phylum = Phylum,	class = Class,	
           order = Order,	family = Family,	genus = Genus,	species = Species)
  # export the file
  write_tsv(taxa_to_FMBN, paste(study$Seq_accn, "_taxa_to_FMBN.txt", sep = ""))
} else {
  cat("no taxa to add to FMBN")
}


# remove some genera ------------------------------------------------------
if(what_data == "B"){
  taxa_ch_lineage <- taxa_ch_lineage |> 
    dplyr::filter(!(str_detect(Genus, "Clostridium|Selenomonas|Incertae_Sedis"))) 
}

if (nrow(taxa_ch_lineage)>0){
  # export the file
  write_tsv(taxa_ch_lineage, paste(study$Seq_accn, "_taxa_ch_lineage.txt", sep = ""))
} else {
  cat(red("no taxa to change"))
} 

# now complete the taxa file with taxonId

taxa_2 <- semi_join(select(taxa_in_FMBN, -id_L6), taxa)
if (dim(taxa_not_in_FMBN)[1]>0){
  taxa_2 <- bind_rows(select(taxa_to_FMBN, taxonId, label), taxa_2)
}

taxa_2 <- full_join(taxa, taxa_2)
# check
beep(which_sound)
cat(red("\nWARNING:\nthis condition must be equal to 0\n"))
anyDuplicated(taxa_2$label)


# fix the edges -----------------------------------------------------------
beep(which_sound)
cat(red("\nWARNING\nBEFORE YOU BEGIN: have you updated the taxa table?\n"))


# edges, may not be needed
# check that all weight sums are 100, must be true
beep(which_sound)
cat(red("\nWARNING\nThis condition must be true\n"))
# must be true 
all(near(pull(edges |> group_by(Run_s) |> summarise(wsum = sum(weight))), 100))

if(dupli_edges && nrow(edges)>nrow(distinct(edges[1:2]))){
  edges <-edges |> 
    group_by(Run_s, s_label) |>
    dplyr::summarise(weight = sum(weight)) |>
    ungroup()
}

# add sampleId, this needs to checked
edges <- left_join(edges, 
                   select(samples, sampleId, Run_s = SRA_run, Source = label_2)) 

# Source is needed for double checking

# add taxonId

edges <- left_join(edges, 
                   select(taxa_2, s_label = label, taxonId))

# add columns for cross checking
edges <- left_join(edges, select(samples, sampleId, Source_samples = label_2)) |>
  left_join(select(taxa_2, taxonId, s_label_taxa = label))

beep(which_sound)
cat(red("\nWARNING:\nthese three conditions must be true\n"))
all(edges$s_label == edges$s_label_taxa)
all(edges$Source == edges$Source_samples)
# check that all weight sums are 100, must be true
all(near(pull(edges |> group_by(Run_s) |> summarise(wsum = sum(weight))), 100))

beep(which_sound)
cat(red("\nWARNING:\nthis condition must be false\n"))
# cross check, must be FALSE
anyNA(edges)

# export file
file_suffix <- if_else(what_data == "B", "_edges_B_to_FMBN.txt", "_edges_F_to_FMBN.txt")
write_tsv(select(edges, sampleId, taxonId, weight, Source, otulabel = s_label), 
          paste(study$Seq_accn, file_suffix, sep = ""))


# edge file has become too large to use it in Excel
# the old edge file (must open both files, bacteria and fungi):
beep(which_sound)
cat(red("\nWARNING\nManually pick the existing edges file for BACTERIA in tab delimited format"))
old_edges_filename_B <- file.choose() # manually pick the old edges in tab delimited format
old_edges_B <- read_tsv(old_edges_filename_B)

beep(which_sound)
cat(red("\nWARNING\nManually pick the existing edges file for FUNGI in tab delimited format"))
old_edges_filename_F <- file.choose() # manually pick the old edges in tab delimited format
old_edges_F <- read_tsv(old_edges_filename_F)

# just a few quick checks
max_sample_edges <- max(old_edges_B$sampleId, old_edges_F$sampleId)
if(!in_FMBN){
  cat(red("\nWARNING:\nthis condition must be true\n"))
  max_samples == max_sample_edges # must be true
  }
max_taxa_edges <- max(old_edges_B$taxonId, old_edges_F$taxonId)

beep(which_sound)
cat(red("\nWARNING:\nthis condition must be true\n"))
max_taxa_edges == max(taxa_in_FMBN$taxonId) # must be true 
# (however, if you have restarted the analysis might be false)
if(!(max_taxa_edges == max(taxa_in_FMBN$taxonId))){
  beep(which_sound)
  cat(red("\nWARNING:\ncarefully check the taxa_in_fmbn, edges and existing edges for coherence\n"))
}

beep(which_sound)
cat(red("\nWARNING:\nthis condition must be false\n"))
ifelse(what_data == "B", anyNA(old_edges_B), anyNA(old_edges_F)) # must be false

edge_file_name <- if_else(what_data == "B",
                          "FMBN_5_edges_B",
                          "FMBN_5_edges_F")
# save the old edges
if(what_data == "B"){
  write_tsv(old_edges_B, str_c(edge_file_name, "_old.txt"))
} else {
  write_tsv(old_edges_F, str_c(edge_file_name, "_old.txt"))  
}

new_edges <- select(edges, sampleId, taxonId, weight)
all_edges <- if(what_data == "B"){
  bind_rows(old_edges_B, new_edges)
} else {
  bind_rows(old_edges_F, new_edges)
  }

write_tsv(all_edges, str_c(edge_file_name, "_new.txt"))



# Package citations -------------------------------------------------------

map(cran_packages, citation)

# Credits and copyright ---------------------------------------------------

# Assume that this is overall under MIT licence

# Copyright 2021, 2022, 2024 Eugenio Parente
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

# This work was partly funded (during years 2024 and 2025) by project
# PRIN2022 PNRR NYCDiversity, Prot. P20229JMMH. Funding was provided by Ministero 
# dell'Università e della ricerca and the European Union Next-GenerationEU 
# (PIANO NAZIONALE DI RIPRESA E RESILIENZA (PNRR) – MISSIONE 4 COMPONENTE 2, 
# INVESTIMENTO 1.4 – D.D. 1032 17/06/2022, CN00000022). 
