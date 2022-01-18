#########1#########2#########3#########4#########5#########6#########7#########8
# merge_phyloseq_objs v1_3
# a script for merging a phyloseq object generated using the DADA2 pipeline with
# phyloseq objects created using the ShinyFMBN app.
# Phylogenetic trees, if present, will be ignored
# A tutorial for the DADA2 pipeline is available here:
# https://benjjneb.github.io/dada2/tutorial.html
# Info on FoodMicrobionet and the Shiny FMBN app can be found at:
# http://www.foodmicrobionet.org
# PLEASE NOTE:
# the script assumes that:
# a. you used SILVA v138.1 for taxonomy assignment (but should work with RDP 
#    trainset as well);
# b. the directory of the script is the working directory;
# c. both objects are saved as in the directory of the script and 
#    they are the only data files in that directory;
# d. the name of the phyloseq object file obtained from ShinyFMBN ends with
#    _physeq.Rdata;
# e. your phyloseq object file is saved as .rds
# The script will check if these conditions hold.
# In general, it may be safer to install the phyloseq package before running
# this script. Go to https://joey711.github.io/phyloseq/install.html for details
#########1#########2#########3#########4#########5#########6#########7#########8
# load packages -----------------------------------------------------------

.cran_packages <- c("tidyverse", "reshape2")
.bioc_packages <- c("BiocManager", "phyloseq")

.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  if(!.inst[1]) install.packages("BiocManager")
  if(any(!.inst[2:length(.inst)])) {
    BiocManager::install(!.inst[2:length(.inst)], ask = F)
  }
}

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

set.seed(100)


# set options -------------------------------------------------------------

change_sample_names <- T
# when set to true with the replace the sample names in your dataset with
# those in the sample name variable (if present in your dataset); the
# occurrence of NAs will be checked. Note that the sample name column
# is renamed in the process

# load files -----------------------------------------------------------
# a function for reading the phyloseq files and checking if they are the correct 
# class
read_phseqs <- function(file_list = list.files()){
  if(all(str_detect(file_list, "merge_phyloseq_objs_", negate = T))){
    warning("The script is not in the working directory, can't load files")
    return()
  }
  # .Rdata files
  Rdata_file <- file_list[str_detect(file_list, "\\.Rdata$|\\.RData$")]
  if(length(Rdata_file)!=1) {
    warning("There must be exactly one .Rdata file in the working directory, stopping without loading")
    return()
    }
  # .rds file
  rds_file <- file_list[str_detect(file_list, "\\.rds$|\\.RDS$")]
  if(length(rds_file)!=1) {
    warning("There must be exactly one .rds file in the working directory, can't identify your phyloseq object")
    return()
    }
  # loading and checking the .Rdata object
  load(Rdata_file)
  if(!exists("physeqdata")) {
    warning("The .Rdata file was not created with ShinyFMBN, stopping without loading")
    return()
    }
  if(class(physeqdata)!="phyloseq") {
    warning("physeqdata is not a phyloseq object, stopping without loading")
    }
  # reading and checking the .rds file
  myphseq <- readRDS(rds_file)
  if(class(myphseq)!="phyloseq") {
    warning("your .rds file is not a phyloseq object, stopping without loading")
    }
  phseqs <- list(FMBN = physeqdata, myphseq = myphseq)
  return(phseqs)
}
# extracts the two phyloseq objects
physeq_list <- read_phseqs()
physeqdata <- physeq_list[[1]]
myphseq <- physeq_list[[2]]


# fixing your phyloseq object ---------------------------------------------

samples <- samples <- as(sample_data(myphseq), "data.frame")

# check which variables might match those from the FMBN samples file
# matching variables: the process is far from perfect, because it is
# impossible to anticipate the content of you sample table

matching <- colnames(samples)[which(colnames(samples) %in% sample_variables(physeqdata))]

# managing columns which are likely to be in the sample table

if(!("label" %in% matching)){
  mypattern <- "Sample_name|Sample_Name|Samplename|SampleName|Sample.Name"
  if(length(which(str_detect(colnames(samples), mypattern)))==0){
    samples$label <- sample_names(myphseq)
  } else {
    colnames(samples)[str_detect(colnames(samples), mypattern)] <- "label"
  }
}

if(!("description" %in% matching)){
  mypattern <- "Description"
  if(length(which(str_detect(colnames(samples), regex(mypattern, ignore_case = T))))==0){
    samples$description <- NA_character_
  } else {
    colnames(samples)[str_detect(colnames(samples), mypattern)] <- "description"
  }
}

if(!("biosample" %in% matching)){
  mypattern <- "bio_sample$|bio_Sample$|Bio_Sample$|biosample$|bioSample$|BioSample$"
  if(length(which(str_detect(colnames(samples), mypattern)))==0){
    samples$biosample <- NA_character_
  } else {
    # only the first matching column is changes
    colnames(samples)[str_detect(colnames(samples), mypattern)] <- "biosample"
  }
}

if(!("SRA_sample" %in% matching)){
  mypattern <- "SRA_sample|SRA_Sample|SRAsample|SRASample"
  if(length(which(str_detect(colnames(samples), mypattern)))==0){
    samples$SRA_sample <- NA_character_
  } else {
    colnames(samples)[str_detect(colnames(samples), mypattern)] <- "SRA_sample"
  }
}

if(!("SRA_run" %in% matching)){
  mypattern <- "SRA_run|SRA_Run|SRArun|SRARun|Run|run"
  if(length(which(str_detect(colnames(samples), mypattern)))==0){
    samples$SRA_run <- NA_character_
  } else {
    colnames(samples)[str_detect(colnames(samples), mypattern)] <- "SRA_run"
  }
}

# managing columns which are unlikely to be in the sample file

if(!("studyId" %in% matching)){
  samples$studyId <- 999
}

if(!("n_reads2" %in% matching)){
  samples$n_reads2 <- NA_real_
}

# fixing NA_character_ columns 

NA_columns_names <- setdiff(sample_variables(physeqdata),names(samples))

samples[,NA_columns_names] <- NA_character_

samples <- samples %>% select_at(.vars = sample_variables(physeqdata))

# handle sample names -----------------------------------------------------

if (change_sample_names){
  # extract the OTU table
  OTUtable <- as(otu_table(myphseq), "matrix")
  # check if the two OTU tables have the same orientation
  # and transpose as needed; the final result is that taxa are rows
  if(!myphseq@otu_table@taxa_are_rows){
    OTUtable <- t(OTUtable)
  }
  # by default sample names have the same order, no need to check
  # unless you are cheating
  row.names(samples) <- samples$label
  colnames(OTUtable) <- samples$label
}


# fix the taxonomy --------------------------------------------------------
# extract the taxonomic table from the object and turn it into a matrix
taxtab <- as(tax_table(myphseq), "matrix")

# transform in a tibble, use the rownames as a column, add species column if needed
taxtab2 <- as_tibble(rownames_to_column(as.data.frame(taxtab, 
                                                      stringsAsFactors = F), "ASV"))
if(!"Species" %in% colnames(taxtab)) taxtab2$Species <- NA_character_

# may be in the future add a variable with the Study, just to merge the database of ASVs

taxtab2[is.na(taxtab2)] <- ""
taxtab2 <- taxtab2 %>%
  mutate(s_label = {
    ifelse(Kingdom == "", "Other",
           ifelse(Phylum == "", Kingdom, 
                  ifelse(Class == "", Phylum,
                         ifelse(Order == "", Class, 
                                ifelse(Family == "", Order, 
                                       ifelse(Genus == "", Family, 
                                              ifelse(Species == "", Genus, paste(Genus, Species, sep =" "))))))))
  })

# fixing bad taxa ---------------------------------------------------------

# first fix taxa with ambiguous labels or potential duplicate labels

# fix Incertae sedis (for Ruminococcaceae and Lachnospiraceae)

pos_to_change <- which(taxtab2$Genus == "Incertae Sedis")
if(length(pos_to_change)>0){
  cat("\nfixing Incertae Sedis in genera")
  taxtab2$Genus[pos_to_change] <- paste(taxtab2$Family[pos_to_change], taxtab2$Genus[pos_to_change], sep = " ")
  taxtab2$s_label[pos_to_change] <- taxtab2$Genus[pos_to_change] }
# fix "Unknown Family"
pos_to_change <- which((taxtab2$Family == "Unknown Family") & (taxtab2$Genus == ""))
if(length(pos_to_change)>0){
  cat("\nfixing Unknown family")
  taxtab2$Family[pos_to_change] <- paste(taxtab2$Order[pos_to_change], taxtab2$Family[pos_to_change], sep = " ")
  taxtab2$s_label[pos_to_change] <- taxtab2$Family[pos_to_change] 
}
# fix s_label == uncultured
pos_to_change <- which(taxtab2$s_label == "uncultured" & taxtab2$Class == "uncultured")
if(length(pos_to_change)>0){
  cat("\nfixing uncultured in Class")
  taxtab2$Class[pos_to_change] <- paste("uncultured", taxtab2$Phylum[pos_to_change], sep = " ")
  taxtab2$s_label[pos_to_change] <- taxtab2$Class[pos_to_change] 
}
pos_to_change <- which(taxtab2$s_label == "uncultured" & taxtab2$Order == "uncultured")
if(length(pos_to_change)>0){
  cat("\nfixing uncultured in Order")
  taxtab2$Order[pos_to_change] <- paste("uncultured", taxtab2$Class[pos_to_change], sep = " ")
  taxtab2$s_label[pos_to_change] <- taxtab2$Order[pos_to_change] 
}
pos_to_change <- which(taxtab2$s_label == "uncultured" & taxtab2$Family == "uncultured")
if(length(pos_to_change)>0){
  cat("\nfixing uncultured in Family")
  taxtab2$Family[pos_to_change] <- paste("uncultured", taxtab2$Order[pos_to_change], sep = " ")
  taxtab2$s_label[pos_to_change] <- taxtab2$Family[pos_to_change] 
}
pos_to_change <- which(taxtab2$s_label == "uncultured" & taxtab2$Genus == "uncultured")
if(length(pos_to_change)>0){
  cat("\nfixing uncultured in Genus")
  taxtab2$Genus[pos_to_change] <- paste("uncultured", taxtab2$Family[pos_to_change], sep = " ")
  taxtab2$s_label[pos_to_change] <- taxtab2$Genus[pos_to_change] 
}

# do further fixes using a lookup table
bad_taxa <- read_tsv(file.path("support","SILVA_bad-taxa.txt"))

if(any(pull(taxtab2, s_label) %in% pull(bad_taxa, s_label))){
  cat("\nfixing bad taxa using lookup table")
  colnames(bad_taxa) <- c("s_label", "new_kingdom", "new_phylum", "new_class", 
                          "new_order", "new_family", "new_genus", "new_species")
  taxtab2 <- left_join(taxtab2, bad_taxa)
  taxtab2 <- taxtab2 %>%
    mutate(Genus = ifelse(is.na(new_genus), Genus, new_genus),
           Family = ifelse(is.na(new_family), Family, new_family),
           Order = ifelse(is.na(new_order), Order, new_order),
           Class = ifelse(is.na(new_class), Class, new_class),
    ) %>%
    select(-(new_kingdom:new_species))
}

# make minor tweaks to handle special cases

if("mitochondria" %in% taxtab2$s_label){
  to_change <- which(taxtab2$s_label == "mitochondria")
  taxtab2$s_label[to_change] <- "Mitochondria"
  print("mitochondria changed")
}

if("Clostridium_sensu_stricto 1" %in% taxtab2$Genus){
  to_change <- which(taxtab2$Genus == "Clostridium_sensu_stricto 1")
  taxtab2$Genus[to_change] <- "Clostridium_sensu_stricto_1"
  taxtab2[to_change] %>%
    mutate(s_label = if_else(s_label == "Clostridium_sensu_stricto 1", 
                           "Clostridium_sensu_stricto_1",
                           s_label))
  print("Clostridium_sensu_stricto 1 changed")
}

if("Clostridium_sensu_stricto 12" %in% taxtab2$Genus){
  to_change <- which(taxtab2$Genus == "Clostridium_sensu_stricto 12")
  taxtab2$Genus[to_change] <- "Clostridium_sensu_stricto_12"
  taxtab2[to_change] %>%
    mutate(s_label = if_else(s_label == "Clostridium_sensu_stricto 12", 
                           "Clostridium_sensu_stricto_12",
                           s_label))
  print("Clostridium_sensu_stricto 12 changed")
}

if("Escherichia-Shigella" %in% taxtab2$Genus){
  to_change <- which(taxtab2$Genus == "Escherichia-Shigella")
  taxtab2$Genus[to_change] <- "Escherichia/Shigella"
  taxtab2[to_change] %>%
    mutate(s_label = if_else(s_label == "Escherichia-Shigella", 
                             "Escherichia/Shigella",
                             s_label))
  print("Escherichia-Shigella changed")
}


# handle changes in species names using a lookup table

lookup <- read_tsv(file.path("support","species_lookup.txt"))

to_change <- which(lookup$id %in% taxtab2$s_label)

taxtab2 <- left_join(taxtab2, lookup, by=c("s_label" = "id"), keep = T)
taxtab2 <- taxtab2 %>%
  mutate(s_label = ifelse(is.na(new_id), s_label, new_id),
         Species = ifelse(is.na(new_species), Species, new_species),
         Genus = ifelse(is.na(new_genus), Genus, new_genus),
         Family = ifelse(is.na(new_family), Family, new_family)
  ) %>%
  select(-(new_id:new_family))

# report changes
if(length(to_change)>0){
  print("The following taxa were changed to make them coherent with FMBN")
  lookup[to_change,]
}

# check the level of aggregation of taxtab2 in the FMBN phyloseq object
taxa_agg <- case_when(
  all(is.na((as(tax_table(physeqdata), "matrix"))[,"order"])) ~ "class",
  all(is.na((as(tax_table(physeqdata), "matrix"))[,"genus"])) ~ "family",
  all(is.na((as(tax_table(physeqdata), "matrix"))[,"species"])) ~ "genus",
  any(!is.na((as(tax_table(physeqdata), "matrix"))[,"species"])) ~ "species"
)

# change variable names to lowercase for compatibility
colnames(taxtab2)[2:8] <- tolower(colnames(taxtab2)[2:8])
# replace "" with NA
replace_empty <- function(x) ifelse(x=="", NA,x)
taxtab2 <- taxtab2 %>% mutate_at(vars(kingdom:species), replace_empty)


taxa_redux <- switch(taxa_agg,
                     "species" = {taxtab2 %>% 
                         mutate(label = s_label) %>%
                         mutate(species = ifelse(is.na(species), NA, s_label))
                         },
                     "genus" = {taxtab2 %>% 
                         mutate(species = NA) %>%
                         mutate(label = ifelse(is.na(genus), s_label, genus))
                       },
                     "family" = {taxtab2 %>% 
                         mutate(genus = NA, species = NA) %>%
                         mutate(label = ifelse(is.na(family), s_label, family))
                       },
                     "class" = {taxtab2 %>% 
                         mutate(order = NA, family = NA, genus = NA, species = NA) %>%
                         mutate(label = ifelse(is.na(class), s_label, class))
                       }
)

# match ASVs and taxonomy
# merge taxonomy into ASV table
seqtab2 <- as_tibble(rownames_to_column(as.data.frame(OTUtable), "ASV"))
seqtab2 <- dplyr::full_join(select(taxa_redux, ASV:species, label), seqtab2) %>%
  select(-ASV) 
# create a long table with "species" labels
seqtab_l <- melt(seqtab2)

# cast and sum
seqtab_ag <- dcast(seqtab_l, label ~ variable, sum)

my_OTU_table <- column_to_rownames(seqtab_ag, var = "label")

# taxonomic table
unique_tax <- taxa_redux %>% select(-ASV, -s_label) %>% distinct()
unique_tax <- unique_tax %>% arrange(label) %>% column_to_rownames("label") 
colnames(unique_tax)[1] <- "domain"
my_tax_table <- as(unique_tax, "matrix")
# build the modified phyloseq object
myphseq_2 <- phyloseq(tax_table(my_tax_table), sample_data(samples),
                      otu_table(my_OTU_table, taxa_are_rows = TRUE))


# create a combined phyloseq object ---------------------------------------

merged_phyloseq <- merge_phyloseq(physeqdata, myphseq_2)


# save the object ---------------------------------------------------------

save(merged_phyloseq, file =file.path("output","merged_phyloseq.Rdata"))

# finis: pro bono, malum


#########1#########2#########3#########4#########5#########6#########7#########8
# Credits and citation ----------------------------------------------------
#########1#########2#########3#########4#########5#########6#########7#########8
# Script created by Eugenio Parente, 2022.                                     #
# The code for installing bioconductor and CRAN packages is taken from:        #
# is derived from https://tinyurl.com/y2pcfdht                                 #
# assume this is under GNU general public licence                              #
# http://www.gnu.org/licenses/gpl-3.0.en.html                                  #
# if you decide to use this script in your work, please mention                #
# the FoodMicrobionet web site                                                 #
# http://www.foodmicrobionet.org                                               #
# and cite relevant publications and datasets therein as appropriate           #
# References for packages used in this script:                                 #
all_packages <- c("base", .cran_packages, .bioc_packages)
map(all_packages, citation)
#########1#########2#########3#########4#########5#########6#########7#########8
# this script was successfully tested on the following systems:

# iMac 21.5 inches late 2013, 2.7 GHz Intel Core i5, quad core, MacOS 10.13.6, R 4.1.2, RStudio 2021.09.0+351 "Ghost Orchid" Release (077589bcad3467ae79f318afe8641a1899a51606, 2021-09-20) for macOS, Safari 15.2  

# MacBook Pro Retina 13\'' 2015, 2.7 GHz Intel Core i5, dual core, MacOS 10.15.7, R 4.1.2, RStudio 2021.09.0+351 "Ghost Orchid" Release (077589bcad3467ae79f318afe8641a1899a51606, 2021-09-20) for macOS, Safari 15.2 

#########1#########2#########3#########4#########5#########6#########7#########8

