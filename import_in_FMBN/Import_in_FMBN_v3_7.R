# Import_in_FMBN ----------------------------------------------------------

# a script for assisting the inclusion of new studies in FMBN# version 3.7 16/11/2021

library(readxl)
library(tidyverse)
library(magrittr)

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

# Studies -----------------------------------------------------------------

# import the study file

# this is bioproject PRJNA643290

study_file <- file_list[!is.na(str_match(file_list, "_study.txt"))]

# open, reorder the columns and save

study <- read_tsv(study_file)
if("tax_database" %in% colnames(study)) {
  study <- select(study, - tax_database)}

#cross-check number of samples, expression must be TRUE
study$samples == n_samples_edge

# ad hoc changed, will not be necessary for more recent versions of the study
# file (after December 2021)

study$primer_f <- "26F4a" # AGAGTTTGATCMTGGCTCAG
study$primer_r <- "534R4" # GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTMTTACCGCGGCNGCTGGCAC contains adapters
study$overlapping <- F
study$paired_end <- T

# make pipeline fields

pip_fields <- tibble(bioinf_software = "R dada2",
                     OTU_picking = "ASV with DADA2",
                     assign_tax_method = "assignTaxonomy()",
                     tax_database = "SILVA v138.1")
empty_fields <- tibble(study = NA, 
                       studyId = NA, 
                       FMBN_version = NA,
                       Seq_accn_link = NA,
                       bioproject = "PRJNA643290",
                       food_group = NA,
                       short_descr = NA,
                       DOI_link = NA,
                       ref_short = NA,
                       year = NA,
                       ref_complete = NA,
                       corr_author_surname = NA, 
                       corr_author_mail = NA)
                       

study <- bind_cols(study, pip_fields, empty_fields) %>%
  select(study,	studyId, FMBN_version, target, region, platform,	
         read_length_bp, seq_center,	bioinf_software,	OTU_picking,	
         assign_tax_method,	tax_database,	Seq_accn,	Seq_accn_link,	bioproject,	
         samples,	food_group, short_descr, DOI_link, ref_short, year, 
         ref_complete, corr_author_surname, corr_author_mail, geoloc, 
         primer_f, primer_r, overlapping, paired_end)

# export the file
write_tsv(study, paste(study$Seq_accn, "_study_FMBN.txt", sep = ""))
         

# samples -----------------------------------------------------------------
# from Excel files FMBNxxx_studies, FMBNxxx_samples
studyId <- "ST172"
max_samples <- 9736
first_sample <- max_samples +1

samples_file <- file_list[!is.na(str_match(file_list, "_samples.txt"))]
samples <- read_tsv(samples_file)



n_samples <- nrow(samples)
rep_NA <- rep(NA, n_samples)
# double check samples number
study$samples == n_samples # must be T

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
                              spoilage = rep_NA
                              )

# this needs to be adjusted manually, to double check inconsistencies
samples <- bind_cols(samples, missing_col_samples) %>%
  mutate(n_reads = n_reads2) %>%
  select(studyId,	sampleId,	label_1,	label_2 = label2,	label_3,	llabel,	s_type,	
         n_reads, n_reads2,	n_issues, foodId,	description,	L1,	L4,	L6,	nature,	
         process, spoilage,	target1,	target2,	biosample,	
         SRA_Sample,	SRA_run, geo_loc_country, geo_loc_continent, lat_lon)

# export the file

write_tsv(samples, paste(study$Seq_accn, "_samples_FMBN.txt", sep = ""))


# taxa --------------------------------------------------------------------

taxa_in_FMBN <- read_excel("taxa_in_fmbn.xlsx")

max_taxaId <- max(taxa_in_FMBN$taxonId)
# check that there are no duplicates
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

nrow(taxa) == n_taxa_edge

# if it fails at first check, need to recheck the original edge file
# still fails, must see which is missing in the edges
length(unique(edges$s_label))
which(duplicated(taxa$label))
# duplication is generally due to "Incertae Sedis" in genus and label
if(length(which(taxa$Genus == "Incertae Sedis"))>0){
  pos_to_change <- which(taxa$Genus == "Incertae Sedis")
  taxa[pos_to_change, ] <- taxa %>% dplyr::filter(Genus == "Incertae Sedis") %>%
    mutate(label = str_c(Family, Genus, sep = " ")) %>%
    mutate(Genus = str_c(Family, Genus, sep = " ")) %>%
    mutate(taxonomy = str_c("k__", Kingdom, "; p__", Phylum, "; c__", Class,
                            "; o__; f__", Family, "; g__", Genus,
                            "; s__")) %>%
    mutate(taxonomy_L6 = str_c("k__", Kingdom, "; p__", Phylum, "; c__", Class,
                            "; o__", Order, "; f__", Family, "; g__", Genus,
                            "; s__")) %>%
    mutate(id = str_c("Root;k__", Kingdom, "__", Phylum, ";c__", Class,
                            ";f__", Family, ";g__", Genus,
                            ";s__")) %>%
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

check_taxa <- check_taxa %>% 
  mutate(issue= (edges != taxa))
all(sort(unique(edges$s_label)) == sort(taxa$label))

# change selected taxa and edges ------------------------------------------

# fix a few taxa (need to be updated every time I spot a new one)
# preserve original versions of taxa and edges
taxa_original <- taxa
edges_original <- edges
# check that all weight sums are 100, must be true
all(near(pull(edges_original %>% group_by(Run_s) %>% summarise(wsum = sum(weight))), 100))

if("mitochondria" %in% taxa$label){
  to_change <- which(taxa$label == "mitochondria")
  taxa$label[to_change,] <- "Mitochondria"
  to_change <- which(edges$s_label == "mitochondria")
  edges$s_label[to_change] <- "Mitochondria"
  print("mitochondria changed")
}

if("Clostridium_sensu_stricto 1" %in% taxa$Genus){
  to_change <- which(taxa$Genus == "Clostridium_sensu_stricto 1")
  taxa$Genus[to_change] <- "Clostridium_sensu_stricto_1"
  taxa[to_change] %>%
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
  taxa[to_change] %>%
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
taxa <- taxa %>%
  mutate(label = ifelse(is.na(new_id), label, new_id),
         Species = ifelse(is.na(new_species), Species, new_species),
         Genus = ifelse(is.na(new_genus), Genus, new_genus),
         Family = ifelse(is.na(new_family), Family, new_family)
  ) %>%
  select(-(new_id:new_family))

edges <- left_join(edges, select(lookup_table, label, new_id), 
                   by=c("s_label" = "label"))
edges <- edges %>%
  mutate(s_label = ifelse(is.na(new_id), s_label, new_id)) %>%
  select(-new_id)

# report changes
if(length(to_change)>0){
  cat("The following taxa were changed to make them coherent with FMBN","\n")
  lookup_table$label[to_change]
} else {
  cat("No changes were made to taxa using the lookup table")
}

# taxa not in FMBN --------------------------------------------------------

taxa_not_in_FMBN <- anti_join(taxa, select(taxa_in_FMBN, label))

# taxa with changed lineages in SILVA 138
taxa_ch_lineage <- inner_join(select(taxa, label:Species, id_L6), select(taxa_in_FMBN, label, id_L6_old = id_L6)) %>%
  dplyr::filter(id_L6 != id_L6_old) %>%
  dplyr::filter(!(label %in% taxa_not_in_FMBN$label))

# double check
in_both <- intersect(taxa_in_FMBN$label, taxa$label)
in_taxa <- setdiff(taxa$label, taxa_in_FMBN$label)
length(in_both) + length(in_taxa) == nrow(taxa)
# must be true 
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
  taxa_to_FMBN <- taxa_not_in_FMBN %>%
    mutate(taxonId = seq(from = max_taxaId+1, 
                         to = max_taxaId+nrow(taxa_not_in_FMBN))) %>%
    select(taxonId,	label,	domain = Kingdom,	phylum = Phylum,	class = Class,	
           order = Order,	family = Family,	genus = Genus,	species = Species)
  # export the file
  write_tsv(taxa_to_FMBN, paste(study$Seq_accn, "_taxa_to_FMBN.txt", sep = ""))
} else {
  cat("no taxa to add to FMBN")
}


# remove some genera ------------------------------------------------------

taxa_ch_lineage <- taxa_ch_lineage %>% 
  dplyr::filter(!(str_detect(Genus, "Clostridium|Selenomonas|Incertae_Sedis"))) 


if (nrow(taxa_ch_lineage)>0){
  # export the file
  write_tsv(taxa_ch_lineage, paste(study$Seq_accn, "_taxa_ch_lineage.txt", sep = ""))
} else {
  cat("no taxa to change")
} 

# now complete the taxa file with taxonId

taxa_2 <- semi_join(select(taxa_in_FMBN, -id_L6), taxa)
if (dim(taxa_not_in_FMBN)[1]>0){
  taxa_2 <- bind_rows(select(taxa_to_FMBN, taxonId, label), taxa_2)
}

taxa_2 <- full_join(taxa, taxa_2)
# check
anyDuplicated(taxa_2$label)

# edges, may not be needed
# check that all weight sums are 100, must be true
all(near(pull(edges %>% group_by(Run_s) %>% summarise(wsum = sum(weight))), 100))

if(dupli_edges && nrow(edges)>nrow(distinct(edges[1:2]))){
  edges <-edges %>% 
    group_by(Run_s, s_label) %>%
    dplyr::summarise(weight = sum(weight)) %>%
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
edges <- left_join(edges, select(samples, sampleId, Source_samples = label_2)) %>%
  left_join(select(taxa_2, taxonId, s_label_taxa = label))
# must be TRUE
all(edges$s_label == edges$s_label_taxa)
all(edges$Source == edges$Source_samples)
# check that all weight sums are 100, must be true
all(near(pull(edges %>% group_by(Run_s) %>% summarise(wsum = sum(weight))), 100))

# cross check, must be FALSE
anyNA(edges)

# export file
write_tsv(select(edges, sampleId, taxonId, weight, Source, otulabel = s_label), 
          paste(study$Seq_accn, "_edges_to_FMBN.txt", sep = ""))

