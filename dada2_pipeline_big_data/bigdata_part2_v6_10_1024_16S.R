################################################################################
# DADA2/Bioconductor pipeline for big data, modified
# part 3, remove chimera, assign taxonomy, save a phyloseq object, prepare files for FMBN
#
# bigdata_part2_v6_9_0824
################################################################################

# This script is designed to process large studies using the
# DADA2 pipeline https://benjjneb.github.io/dada2/tutorial.html with options
# for large studies https://benjjneb.github.io/dada2/bigdata.html

# This version of the script (v6_9) combines sequence tables and performs 
# taxonomic assignment and post processing, following part 1

# to use this script follow the instructions in script bigdata_part1 and run  
# getHeader and bigdata_part1 (as many times as there are groups in the identifier file)

# You can easily adapt the script to your own data

# Run the script one section at a time

# load packages -----------------------------------------------------------

.cran_packages <- c("tidyverse", "parallel", "gridExtra", "knitr", "stringr", 
                    "reshape2", "beepr", "tictoc", "logr")
.bioc_packages <- c("BiocManager","dada2", "phyloseq", "BiocStyle", 
                    "DECIPHER", "phangorn")

.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  if(!.inst[1]) {
    install.packages("BiocManager")
    .inst <- .bioc_packages %in% installed.packages()
  }
  if(any(!.inst[2:length(.inst)])) {
    BiocManager::install(.bioc_packages[!.inst], ask = F)
  }
}

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

sessionInfo()

# other setup operations ------------------------------------------------

opar <- par(no.readonly=TRUE) 
par(ask=F) 
# set.seed(1234) 
# plays audio notifications if TRUE
play_audio <- T
sound_n = 6 # an integer from 1 to 11, sets the notification sound
# keeps record of duration of important steps if TRUE
keep_time <- T
# verbose output: will print additional objects and messages if TRUE
verbose_output <- T
# do you want to use a log to store the options you used in processing?
# the default is T with verbose_output <- T otherwise you have to set it
# manually
use_logr <- ifelse(verbose_output, T, F)
# use_logr <-T use this line to change it manually

if(play_audio) beep(sound = sound_n) # the notification you will hear
# the following command detects the number of cores on UNIX/MacOS

nc <- parallel::detectCores(logical = F) # to detect physical cores in MacOS
set.seed(100)


# info for file saving

# load data from previous session
gen_data <- readRDS(file = "gen_data.RDS")

Study <- gen_data$Study
target <- gen_data$target
region <- gen_data$region
seq_accn <- Study
DOI <- gen_data$DOI

# information on the platform and arrangement
data_type <- gen_data$data_type
platform <- gen_data$platform 
overlapping <- gen_data$overlapping
target1 <- gen_data$target1
target2 <- region
paired_end <- gen_data$paired_end 
primer_f <- gen_data$primer_f 
primer_r <- gen_data$primer_r 
FWD <- gen_data$primerf_seq
REV <- gen_data$primerr_seq

# information on merging
merge_option <- gen_data$merge_option


# merge multiple runs -----------------------------------------------------
# note for self: maybe use a functional instead
# get the sequence files
seq_files <- list.files()[grepl("seqtab_f_", list.files())]
seq_tables <- str_remove(seq_files, ".Rdata")
for(i in seq_along(seq_files)) {
  load(seq_files[i])
  assign(seq_tables[i],seqtab_f)
  }

# pass the names to the mergeSequenceTables function

st.all <- mergeSequenceTables(tables = lapply(seq_tables, get))


# get track files
track_files <- list.files()[grepl("track_", list.files())]

track_list <-lapply(track_files, readRDS)
track <- base::do.call("rbind", track_list)

write_tsv(as.data.frame(track), str_c("track_all_",Study,".txt"))


# remove bimera -----------------------------------------------------------

# may be very slow for large studies, you may want to check with a shorter version
# seqtab.nochim <- removeBimeraDenovo(seqtab_f[,1:10000])
dim(st.all)
seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
beep(sound=6)
# better save workspace here (and then remove large objects to recover memory)

save.image(file = str_c(Study,".Rdata"))
rm(list = ls(pattern = "seqtab_f_"))

# recheck length distribution and see if singletons are still there
# look at the abundance distribution
data.frame(nseqs = colSums(seqtab.nochim)) |>
  ggplot() +
  geom_histogram(mapping = aes(x= log10(nseqs))) +
  labs(
    title = "sequence length distribution", 
    x = "log10(seq_count)", 
    y = "number of individual sequences"
  )

# which ASVs are singletons or doubletons?
single_double <- which(colSums(seqtab.nochim)<=2)
length(single_double)/ncol(seqtab.nochim)
# fraction of singletons+doubletons here 5% (remember, these are ASVs, not OTUs) 
singletons <- which(colSums(seqtab.nochim)<=1)
length(singletons)
length(singletons)/ncol(seqtab.nochim)
doubletons <- which(colSums(seqtab.nochim)==2)
length(doubletons)

# If you want to remove singletons and doubletons set remove_s_d to T
remove_s_d <- F
seqtab.nochim.all <- seqtab.nochim
if(remove_s_d){
  seqtab.nochim <- seqtab.nochim[,-single_double]
}

dim(seqtab.nochim)
# check which is the abundance of chimeras removed
sum(seqtab.nochim)/sum(st.all)
# check the distribution of removals
summary(rowSums(seqtab.nochim)/rowSums(st.all))

# save and remove unneeded objects 
save.image(file = str_c(Study,".Rdata"))
rm(seqtab.nochim.all, st.all)

# complete the track file -------------------------------------------------

track <- cbind(track, rowSums(seqtab.nochim))
colnames(track)[6] <- "nonchim"
track2 <- as.data.frame(track) %>%
  rownames_to_column(var = "sample_name") %>%
  dplyr::mutate(
    prop_filt = filtered/input,
    prop_denoised = denoised/input,
    prop_merged = merged/input,
    prop_tabled = tabled/input,
    prop_nonchim = nonchim/input
  )

# flagging the samples with excessive loss of features or too few features
track2 <- track2 %>% dplyr::mutate(
  high_chim = dplyr::if_else((tabled-nonchim)/tabled >0.30, 1, 0),
  high_filt_loss = dplyr::if_else(prop_filt < 0.55, 1, 0),
  low_feat = dplyr::if_else(prop_nonchim < 0.10, 1, 0),
  low_seq_n = dplyr::if_else(nonchim < 5000, 1, 0)
)
track2$n_issues <- rowSums(dplyr::select(track2, high_chim:low_seq_n))

write_tsv(as.data.frame(track2), str_c("track_all_",Study,".txt"))


# assign taxonomy ---------------------------------------------------------

# set the directory for taxonomy databases
taxdb_dir <- file.path("..","tax_db") # change this if the tax databases are elsewhere
list.files(taxdb_dir)

# assignment with SILVA
RC <- F # option for trying reverse complement, false by default

# Loading files needed to manage taxonomy, note locations
ref_fasta <- file.path(taxdb_dir, "silva_nr99_v138_1_train_set.fa")
sp_ass_SILVA <- file.path(taxdb_dir, "silva_species_assignment_v138_1.fa")


# defining the function for assigning taxonomy ----------------------------

# function for assigning taxonomy and returning a list of objects
# needed to avoid copy and paste in the loop

Assign_taxonomy <- function(seqtab, paired_end = T, overlapping = T, rc = RC){
  cat("Assigning taxonomy at the genus level with SILVA v138...", "\n")
  taxtab <- assignTaxonomy(seqtab, refFasta = ref_fasta, multithread = TRUE, tryRC = rc)
  cat("done...", "\n")
  # do species assignment (only if overlapping or not paired end)
  if(!paired_end | overlapping){
    cat("Assigning taxonomy at the species level with SILVA v138...", "\n")
    taxtab <- addSpecies(taxtab, sp_ass_SILVA, tryRC = rc)
    cat("done...", "\n")
  }
  taxtab2 <- as_tibble(rownames_to_column(as.data.frame(taxtab, 
                                                        stringsAsFactors = F), "ASV"))
  if(!"Species" %in% colnames(taxtab)) taxtab2$Species <- NA_character_
  
  # maybe in the future add a variable with the Study, just to merge the database of ASVs
  
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
  
  # build the list to return
  cat("building the list to return", "\n")
  tax_tab_list = list(
    taxtab_slot = taxtab,
    taxtab2_slot = taxtab2
  )
  cat("done...", "\n")
  return(tax_tab_list)
}

split_point <- 2500
# if <2500 ASVs, run in one set, else split and use a loop
# else run only once
# roughly for short sequences, like V4 you can set to 5000, while for V1-V3 
# and V3-V4 2500 is safer with 8 Gb RAM

split <- F
if(dim(seqtab.nochim)[2]>=split_point) split <-T


# set options for taxonomy assignment -------------------------------------

# check if this matches your dataset
pend <- paired_end # check if T or F
ovlp <- overlapping

if(!split){
  cat("assigning taxonomy...", "\n")
  
  taxtab_list <- Assign_taxonomy(seqtab = seqtab.nochim, 
                                 paired_end = pend, overlapping = ovlp)
  if(mean(is.na(taxtab_list$taxtab_slot[,2]))>0.2){
    RC<<-T
    taxtab_list <- Assign_taxonomy(seqtab = seqtab_nochim, 
                                   paired_end = pend, overlapping = ovlp, rc = RC)
  }
  cat("done...","\n")
} else {
  # this is not tested yet
  maxcolseqtab <- dim(seqtab.nochim)[2]
  first_column <- 1
  iterations <- ifelse((maxcolseqtab-first_column)%%split_point,
                       (maxcolseqtab-first_column)%/%split_point+1,
                       (maxcolseqtab-first_column)%/%split_point)
  # iterate over the sub-groups 
  for(i in 1:iterations){
    cat(str_c("iteration", i, "of", iterations, sep = " "),"\n")
    from_c <- ifelse(i==1, first_column, first_column+split_point*(i-1))
    to_c <- ifelse((from_c+split_point)<=maxcolseqtab, (from_c+split_point-1), maxcolseqtab)
    cat(str_c("ASV", from_c, "to", to_c, sep = " "),"\n")
    seqtab_nochim_temp <- seqtab.nochim[,from_c:to_c]
    cat("Assigning taxonomy...","\n")
    taxtab_list_temp <- Assign_taxonomy(seqtab = seqtab_nochim_temp, 
                                        paired_end = pend, overlapping = ovlp, rc = RC)
    if(i == 1 && mean(is.na(taxtab_list_temp$taxtab_slot[,2]))>0.2){
        RC<<-T
        taxtab_list_temp <- Assign_taxonomy(seqtab = seqtab_nochim_temp, 
                                            paired_end = pend, overlapping = ovlp, rc = RC)
    }
    
    cat("done...","\n")
    if(i==1) {
      taxtab_list <- taxtab_list_temp

    } else {
      taxtab_list$taxtab_slot <- rbind(taxtab_list$taxtab_slot, taxtab_list_temp$taxtab_slot)
      taxtab_list$taxtab2_slot <- rbind(taxtab_list$taxtab2_slot, taxtab_list_temp$taxtab2_slot)
      gc() # do garbage collection
    }
    # save the temporary list
    cat("saving the list, iteration ", as.character(i), "\n")
    saveRDS(taxtab_list, str_c("taxtablist_", as.character(i), ".RDS"))
    cat("done...","\n")
  }
  
# reclaim memory

gc() # do garbage collection

# better save here

save.image(file = str_c(Study,".Rdata"))
beep(sound=6)  
}

# fixing the output

taxtab <- taxtab_list$taxtab_slot
taxtab2 <- taxtab_list$taxtab2_slot

if(!split) {
  rm(taxtab_list, taxtab_list_temp, Assign_taxonomy)
} else {
  rm(taxtab_list, Assign_taxonomy) 
}
gc()

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
bad_taxa <- read_tsv(file.path(taxdb_dir,"SILVA_bad-taxa.txt"))

if(any(pull(taxtab2, s_label) %in% pull(bad_taxa, s_label))){
  cat("\nfixing bad taxa using lookup table\n")
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



# save a smaller version of workspace which will be needed later
save(taxtab, seqtab.nochim, track2, Study, target, region, 
     seq_accn, DOI, primer_f, primer_r, target1, target2, taxtab2, 
     file = str_c(Study,"_small.Rdata"))



# Build phylogenetic tree ---------------------------------------------

# decipher can align only ~46k seqs
# https://www.bioconductor.org/packages/3.7/bioc/vignettes/DECIPHER/inst/doc/ArtOfAlignmentInR.pdf
dim(seqtab.nochim)[2]
# if this is >46k will crash
# but it is already prohibitive with 20k and usually not worth the effort

dotree <- F # avoid if overlapping == F

if (dotree) {
  seqs <-
    dada2::getSequences(seqtab.nochim) # the collapse option is very interesting
  names(seqs) <- seqs # This propagates to the tip labels of the tree
  # run alignment from the decipher package 11:47 end 11:48
  alignment <-
    DECIPHER::AlignSeqs(DNAStringSet(seqs),
                        anchor = NA,
                        processors = nc)
  
  # The phangorn R package is then used to construct a phylogenetic tree.
  # Here we first construct a neighbor-joining tree, and then fit a GTR+G+I
  # (Generalized time-reversible with Gamma rate variation)
  # maximum likelihood tree using the neighbor-joining tree as a starting point.
  # tansform in phydat object
  phang.align <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
  # create distance matrix
  cat("Creating distance matrix...", "\n")
  dm <- phangorn::dist.ml(phang.align)
  cat("Performing Neighbor joining...", "\n")
  # perform Neighbor joining
  treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order
  cat("Calculating internal maximum likelihood...", "\n")
  # internal maximum likelihood for tree
  fit = phangorn::pml(treeNJ, data = phang.align)
  Sys.sleep(5)
  fitGTR <- update(fit, k = 4, inv = 0.2)
  Sys.sleep(5)
  # this is the step taking the longest time
  cat("Optimization, please be patient (with >1000 seqs you are better off doing this overnight)...","\n")
  fitGTR <- optim.pml(
    fitGTR,
    model = "GTR",
    optInv = TRUE,
    optGamma = TRUE,
    rearrangement = "stochastic",
    control = pml.control(trace = 0)
  )
  detach("package:phangorn", unload = TRUE)
}
beep(sound = 6)

save.image(file = str_c(Study,"_small.Rdata"))

# Combine data into a phyloseq object ------------------------------------- 

# path to the metadata file (needs to be adapted)
metadata_path <- file.path("data", "metadata", "SraRunTable.txt") # "data/metadata/SraRunInfo.txt"
# alternatives are:
# "data/metadata/SraRunInfo.txt" 
# "data/metadata/SraRunTable.txt.csv"
# "data/metadata/SraRunTable.txt"
samdf <- read_tsv(metadata_path) 
if(ncol(samdf)==1) samdf <- read_csv(metadata_path)

# check if any cols in seqtab have 0 sums
seq_sums <- rowSums(seqtab.nochim)
runs_to_keep <- which(seq_sums>0)
runs <- rownames(seqtab.nochim)[runs_to_keep]
if(length(runs_to_keep) < dim(seqtab.nochim)[1]) {
  seqtab.nochim <- seqtab.nochim[runs_to_keep,]
  samdf <- dplyr::filter(samdf, Run %in% runs)
}

# sanity check
all(rownames(seqtab.nochim) %in% samdf$Run) # must be TRUE
rownames(seqtab.nochim) <- sapply(strsplit(rownames(seqtab.nochim), "_", fixed = T), `[`, 1)
all(rownames(seqtab.nochim) %in% samdf$Run) # TRUE
# add final number of sequences and number of issues from track2
# add number of sequences to track2
seqs <- as.data.frame(rowSums(seqtab.nochim))
seqs <- rownames_to_column(seqs, var = "Run")
colnames(seqs)[2] <-"seqs"
track2 <- left_join(track2, select(data = seqs, sample_name = Run, seqs))

samdf <- full_join(samdf, select(track2, Run = sample_name, n_reads2 = data, 
                                 n_issues))
samdf <- as.data.frame(samdf)
rownames(samdf) <- samdf$Run

track2 <- mutate(track2, prop_tax_filt = data/nonchim)

rownames(samdf) <- samdf$Run

# a transposed sequence table
seqtab_t <- t(seqtab.nochim)
dim(seqtab_t)

# may have dropped a few runs because they lacked seqs
samdf_2 <- samdf %>% dplyr::filter(Run %in% colnames(seqtab_t))

# just checking
# colnames(seqtab_t) == row.names(samdf_2)
all(colnames(seqtab_t) == row.names(samdf_2))
# see if needs reordering
# sort(colnames(seqtab_t)) == sort(row.names(samdf_2))
all(sort(colnames(seqtab_t)) == sort(row.names(samdf_2)))

# combine in a phyloseq object
if(exists("fitGTR", mode = "list")){
  myphseq <- phyloseq(tax_table(taxtab), sample_data(samdf_2),
                      otu_table(seqtab_t, taxa_are_rows = TRUE),
                      phy_tree(fitGTR$tree))
} else {
  myphseq <- phyloseq(tax_table(taxtab), sample_data(samdf_2),
                      otu_table(seqtab_t, taxa_are_rows = TRUE))}


# save the object for further analysis
saveRDS(myphseq,str_c("data/", Study, "_ps.rds"))

# save the workspace
save.image(file = str_c(Study,"_small.Rdata"))

# prepare files for FMBN ----------------------------------------------

# The following instructions are needed to prepare objects which, after further
# processing using the Import_in_FMBN_xx.R script, will be used to generate
# .txt files for import i FMBN tables

# samples -----------------------------------------------------------------
# loads the phyloseq object if it does not exist
if (!exists("myphseq")) {myphseq <- readRDS(str_c("data/", Study, "_ps.rds", sep = ""))}

# prep the sample table
# extract the sample data to a data frame
samples <- as(sample_data(myphseq), "data.frame")
# extract the tax table to a matrix (however, I will be using taxtab2)
ttab <- as(tax_table(myphseq), "matrix")

##############################
# check the runInfotable and adapt these statements
##############################
#  study information (from sample table) ---------------------------------------
n_samples <- nrow(samples)
instrument <- unique(samples$Instrument)[1]  # in some cases length >1
seq_center <- unique(samples$Center_Name) 
if(is_null(seq_center)) seq_center <- unique(samples$Center.Name)
if(is_null(seq_center)) seq_center <- unique(samples$CenterName)
# extract the average read length as an integer
read_length <- round(mean(nchar(rownames(ttab)), na.rm = T)) 

# adapt this or set manually a comma delimited string of countries.
loc_list <- ifelse("geo_loc_name_country" %in% colnames(samples),
                   str_flatten(pull(distinct(samdf, geo_loc_name_country)), collapse =", "),
                   NA_character_)

# put together and save study info
study <- tibble(target = target, region = region, platform = instrument,
                read_length_bp = read_length, seq_center = seq_center,
                tax_database = "SILVA v138_1", Seq_accn = seq_accn,
                samples = n_samples, DOI = DOI, geoloc = loc_list, 
                primer_f, primer_r, overlapping = overlapping, 
                paired_end = paired_end)
# saves study info
write_tsv(study, str_c(Study,"_study.txt"))


# format sample information for FMBN -------------------------------------------

# this needs to be adapted for each study
# check naming of the geoloc info

samples <- samples %>%
  mutate(description = str_c(env_local_scale, env_medium, Library.Name, sep =", "))
samples <- samples %>%
  mutate(Sample_Name = Run) 

# information of geoloc (and names of the field) is very inconsistent:
# check the info in your sample metadata and adatp these commands
# use these if part or all of the geolocation information is missing
# samples$geo_loc_name_country <- NA_character_
# samples$geo_loc_name_country_continent <- NA_character_
# samples$lat_lon <- NA_character_


# create label2 (to avoid numbers as first char.; s. can be removed later with
# tidyr::separate)
# needs to be adjusted ad hoc especially for lat_lon
if(data_type == "sra"){
  samples <- samples %>%
    mutate(label2 = str_c("s.",Sample_Name), target1 = target1, 
           target2 = target2) %>%
    select(label2, n_reads2, n_issues, description, target1, target2, 
           biosample = BioSample, SRA_Sample = BioSample, SRA_run = Run, 
           geo_loc_country = geo_loc_name_country, 
           geo_loc_continent = geo_loc_name_country_continent, lat_lon)
} else {
  samples <- samples %>%
    mutate(label2 = str_c("s.",Sample_Name), target1 = target1, 
           target2 = target2, SRA_Sample = NA_character_, 
           SRA_run = NA_character_) %>%
    select(Run, label2, n_reads2, n_issues, description, target1, target2, 
           biosample = Sample_code, geo_loc_country = Country, 
           geo_loc_continent = Continent, lat_lon, Sample_Name)
}


# save the sample information
write_tsv(samples, str_c(Study,"_samples.txt"))


# extract unique taxa -----------------------------------------------------


# a tibble with unique elements for taxonomy
# some changes are necessary for coherence with FMBN taxonomy
unique_tax <- taxtab2 %>% select(-ASV) %>% distinct()


# if the tax database is silva v138 class changes are not needed

if(!str_detect(study$tax_database, "v138")) {
  taxa <- unique_tax %>%
    mutate(sp_label = ifelse(Species!="",str_c(Genus, Species),"")) %>%
    mutate(Class = ifelse(Class == "Acidomicrobia", "Acidomicrobiia", Class)) %>%
    mutate(Class = ifelse(Class == "Coriobacteria", "Coriobacteriia", Class)) %>%
    mutate(Class = ifelse(Class == "Flavobacteria", "Flavobacteriia", Class)) %>%
    mutate(Class = ifelse(Class == "Sphingobacteria", "Sphingobacteriia", Class)) %>%
    mutate(Class = ifelse(Class == "Fusobacteria", "Fusobacteriia", Class)) %>%
    mutate(id = ifelse(Kingdom =="", "Root;k__;c__;f__;g__;s__",
                       str_c("Root;k__", Kingdom, "__", Phylum, ";c__", 
                             Class, ";f__", Family, ";g__", Genus, ";s__", sp_label)),
           id_L6 = ifelse(Kingdom =="", "Root;k__;c__;o__;f__;g__;s__",
                          str_c("Root;k__", Kingdom, "__", Phylum, ";c__", Class, 
                                ";o__", Order,  ";f__", Family, ";g__", Genus, 
                                ";s__", sp_label))) %>%
    mutate(sp_label = ifelse(Species!="",str_c(Genus, Species, sep = " "),"")) %>%
    mutate(taxonomy = str_c("k__", Kingdom, "; p__", Phylum, "; c__", Class,
                            "; o__; f__", Family, "; g__", Genus, "; s__", sp_label),
           taxonomy_L6 = str_c("k__", Kingdom, "; p__", Phylum, "; c__", Class,
                               "; o__", Order,  "; f__", Family, "; g__", Genus, 
                               "; s__", sp_label)) %>%
    mutate(Species = sp_label) %>%
    select(id, label = s_label, Kingdom:Species, taxonomy, id_L6, taxonomy_L6)
  
  # found a bug, need to fix labels for some taxa
  taxa <- taxa %>% mutate(label = case_when(
    Class == "Acidobacteria (class)" & label == "Acidobacteria" ~ "Acidobacteria (class)",
    TRUE ~ label
  ))
  taxtab2 <- taxtab2 %>% mutate(s_label = case_when(
    Class == "Acidobacteria (class)" & s_label == "Acidobacteria" ~ "Acidobacteria (class)",
    TRUE ~ s_label
  ))
} else {
  taxa <- unique_tax %>%
    mutate(Genus = ifelse(Genus == "Incertae_sedis", str_c(Family, Genus, sep = "_"), Genus)) %>%
    mutate(sp_label = ifelse(Species!="",str_c(Genus, Species),"")) %>%
    mutate(id = ifelse(Kingdom =="", "Root;k__;c__;f__;g__;s__",
                       str_c("Root;k__", Kingdom, "__", Phylum, ";c__", 
                             Class, ";f__", Family, ";g__", Genus, ";s__", sp_label)),
           id_L6 = ifelse(Kingdom =="", "Root;k__;c__;o__;f__;g__;s__",
                          str_c("Root;k__", Kingdom, "__", Phylum, ";c__", Class, 
                                ";o__", Order,  ";f__", Family, ";g__", Genus, 
                                ";s__", sp_label))) %>%
    mutate(sp_label = ifelse(Species!="",str_c(Genus, Species, sep = " "),"")) %>%
    mutate(taxonomy = str_c("k__", Kingdom, "; p__", Phylum, "; c__", Class,
                            "; o__; f__", Family, "; g__", Genus, "; s__", sp_label),
           taxonomy_L6 = str_c("k__", Kingdom, "; p__", Phylum, "; c__", Class,
                               "; o__", Order,  "; f__", Family, "; g__", Genus, 
                               "; s__", sp_label)) %>%
    mutate(Species = sp_label) %>%
    select(id, label = s_label, Kingdom:Species, taxonomy, id_L6, taxonomy_L6)
}


write_tsv(taxa, str_c(Study, "taxaFMBN.txt"))         

# prepare OTU and edge tables ---------------------------------------------

# merge taxonomy into sequence table
seqtab2 <- as_tibble(rownames_to_column(as.data.frame(otu_table(myphseq)), "ASV"))
seqtab2 <- dplyr::full_join(select(taxtab2, ASV, s_label), seqtab2) %>%
  select(-ASV) 


# create a long table with "species" labels
seqtab_l <- pivot_longer(seqtab2, cols = 2:ncol(seqtab2), names_to = "variable")

# adds a genus label and the full taxonomy
seqtab_lg <- dplyr::full_join(seqtab_l, unique_tax) %>% 
  mutate(g_label = ifelse(Genus == "", s_label, Genus)) 

#the edge table
edge_table <- seqtab_lg %>% 
  group_by(variable, s_label) %>%
  dplyr::summarise(seq_sums = sum(value)) %>%
  dplyr::rename(Run_s = variable) %>%
  ungroup() %>%
  group_by(Run_s) %>%
  mutate(weight = 100*seq_sums/sum(seq_sums)) %>%
  dplyr::filter(weight>0) %>%
  ungroup()

# get sums at the species level and pivot wider, absolute frequencies

seqtab_ags <- seqtab_lg %>%
  group_by(variable, s_label) %>%
  dplyr::summarise(seq_sums = sum(value)) %>%
  pivot_wider(id_cols = s_label, names_from = variable, values_from = seq_sums)

# same genus level
seqtab_agg <- seqtab_lg %>%
  group_by(variable, g_label) %>%
  dplyr::summarise(seq_sums = sum(value)) %>%
  pivot_wider(id_cols = g_label, names_from = variable, values_from = seq_sums)

# creates tables of proportions
seqtab_agsp <- seqtab_ags %>% mutate_if(is.numeric, function(.)(./sum(., na.rm = T)))
seqtab_aggp <- seqtab_agg %>% mutate_if(is.numeric, function(.)(./sum(., na.rm = T)))

# saving objects ----------------------------------------------------------

write_tsv(seqtab_ags, str_c(Study,"seqtab_ags.txt"))
write_tsv(seqtab_agg, str_c(Study,"seqtab_agg.txt"))
write_tsv(seqtab_agsp, str_c(Study,"seqtab_agsp.txt"))
write_tsv(seqtab_aggp, str_c(Study,"seqtab_aggp.txt"))
write_tsv(edge_table, str_c(Study,"edge_table.txt"))

# save the workspace
save.image(file = str_c(Study,"_small.Rdata"))

# create a list for reprocessing taxonomy if needed

# assembling and saving the list -----------------------------------------------------

# checking if everything which is needed is available ---------------------

check_list <- c(
  exists("Study"),
  exists("study"),
  exists("myphseq"),
  exists("overlapping"),
  exists("paired_end"),
  exists("RC")
)


mylist <- list(Study_accn = Study,
               study_df = study,
               overlap = overlapping,
               pend = paired_end,
               physeq = myphseq,
               rev_compl = RC)


if(all(check_list)){
  cat("\nAll needed objects available and ready to process\n")
  saveRDS(mylist, file = str_c(Study, "_mindata.RDS"))
  cat("\nSaved data for",Study,"\n")
} else {
  cat("\nOne or more of the objects you need is missing, check your data before proceeding\n")
}

# save the workspace
save.image(file = str_c(Study,"_small.Rdata"))


# Package citations -------------------------------------------------------

map(c(.cran_packages, .bioc_packages), citation)



# Credits and copyright ---------------------------------------------------

# Most of the script is taken from https://benjjneb.github.io/dada2/tutorial.html
# or https://benjjneb.github.io/dada2/bigdata.html
# with changes and adaptations

# Assume that this is overall under MIT licence

# Copyright 2021, 2022, 2024, Eugenio Parente
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



