#########1#########2#########3#########4#########5#########6#########7#########8
# make_seq_df.R a script for creating a data frame with ref seqs  --------
#########1#########2#########3#########4#########5#########6#########7#########8

# the script takes as an input fasta files with taxonomic reference databases
# and convert them in a df
# the fasta files must be in a folder named taxdb one level up from the project
# or working directory
# the SILVA v138.1 reference databases used in this script were downloaded from
# https://doi.org/10.5281/zenodo.4587954

library(phylotools)
library(tidyverse)


# load the reference database ---------------------------------------------

taxdb_dir <- "../tax_db" # change this if the tax databases are elsewhere
list.files(taxdb_dir)
ref_database_genus <- "silva_nr99_v138_1_train_set.fa"
ref_database_genus_species <- "silva_nr99_v138_1_wSpecies_train_set.fa"

ref_fasta_SILVA_genus <- file.path(taxdb_dir, ref_database_genus)
ref_fasta_SILVA_genus_species <- file.path(taxdb_dir, ref_database_genus_species)

SILVA_genus_df <- read.fasta(file = ref_fasta_SILVA_genus, clean_name = FALSE)
# a few changes for compatibility with FoodMicrobionet
SILVA_genus_df_2 <- SILVA_genus_df %>%
  separate(seq.name, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
           sep = ";", remove = F) %>%
  mutate(Species = NA_character_) %>%
  mutate(seq.name = {
    ifelse(is.na(Kingdom), "Other",
           ifelse(is.na(Phylum), Kingdom, 
                  ifelse(is.na(Class), Phylum,
                         ifelse(is.na(Order), Class, 
                                ifelse(is.na(Family), Order, 
                                       ifelse(is.na(Genus), Family, 
                                              ifelse(is.na(Species), Genus, paste(Genus, Species, sep =" ")))))))) 
  }) %>%
  dplyr::filter(seq.name != "")
write_tsv(SILVA_genus_df_2, "SILVA_genus.txt")
rm(SILVA_genus_df, SILVA_genus_df_2)
gc()

# same with species references
SILVA_genus_species_df <- read.fasta(file = ref_fasta_SILVA_genus_species, clean_name = FALSE)
SILVA_genus_species_df_2 <- SILVA_genus_species_df %>%
  separate(seq.name, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", remove = F) %>%
  mutate(seq.name = {
    ifelse(is.na(Kingdom), "Other",
           ifelse(is.na(Phylum), Kingdom, 
                  ifelse(is.na(Class), Phylum,
                         ifelse(is.na(Order), Class, 
                                ifelse(is.na(Family), Order, 
                                       ifelse(is.na(Genus), Family, 
                                              ifelse(is.na(Species), Genus, paste(Genus, Species, sep =" ")))))))) 
  }) %>%
  dplyr::filter(seq.name != "")

write_tsv(SILVA_genus_species_df_2, "SILVA_genus_species.txt")
rm(SILVA_genus_species_df, SILVA_genus_species_df_2)
gc()