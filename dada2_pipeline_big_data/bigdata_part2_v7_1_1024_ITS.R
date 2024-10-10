################################################################################
# DADA2/Bioconductor pipeline for big data, modified
# part 3, remove chimera, assign taxonomy, save a phyloseq object, prepare files for FMBN
#
# bigdata_part2_v7_0724_ITS
################################################################################

# This script is designed to process large studies using the
# DADA2 pipeline https://benjjneb.github.io/dada2/tutorial.html with options
# for large studies https://benjjneb.github.io/dada2/bigdata.html

# This version of the script (v7) includes options
# for single end/paired end data sets obtained with Illumina or 454 or
# Ion Torrent, is adapted for taxonomy assignment with UNITE,
# and will perform steps up to the creation of the sequence table
# For each iteration of the script over the identifiers.txt table, created with
# the getHeader script, a new sequence table is created and saved from the
# corresponding group fo sequences. These are then assembled and processed
# using the bigdata_part2 script which carries out further processing needed to ready 
# the output for import into FoodMicrobionet

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

sessionInfo()


# load data from previous session -----------------------------------------

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
if(play_audio) beep(sound=sound_n)

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

if(keep_time) tic("\nassign taxonomy")
# set the directory for taxonomy databases

RC <-F # try reverse complement of asignTaxonomy(), will be changed if needed

# UNITE
taxdb_dir <- file.path("..","tax_db") # change this if the tax databases are elsewhere
ref_fasta <- file.path(taxdb_dir, "sh_general_release_dynamic_s_04.04.2024.fasta")

# defining the function for assigning taxonomy ----------------------------

# function for assigning taxonomy and returning a list of objects
# needed to avoid copy and paste in the loop

Assign_taxonomy_ITS <- function(seqtab, paired_end = T, overlapping = T, rc = RC){
  cat("Assigning taxonomy with UNITE...", "\n")
  taxtab <- assignTaxonomy(seqtab, refFasta = ref_fasta, multithread = TRUE, tryRC = rc)
  cat("done...", "\n")
  taxtab2 <- as_tibble(rownames_to_column(as.data.frame(taxtab, 
                                                        stringsAsFactors = F), "ASV"))
  if(!"Species" %in% colnames(taxtab)) taxtab2$Species <- NA_character_
  
  # clean taxa
  taxtab2 <- taxtab2 %>%
    mutate(Kingdom = str_remove(Kingdom, "k__"),
           Phylum = str_remove(Phylum, "p__"),
           Class = str_remove(Class, "c__"),
           Order = str_remove(Order, "o__"),
           Family = str_remove(Family, "f__"),
           Genus = str_remove(Genus, "g__"),
           Species = str_remove(Species, "s__"))
  
  taxtab2[is.na(taxtab2)] <- ""
  taxtab2 <- taxtab2 %>%
    mutate(s_label = {
      ifelse(Kingdom == "", "Other",
             ifelse(Phylum == "", Kingdom, 
                    ifelse(Class == "", Phylum,
                           ifelse(Order == "", Class, 
                                  ifelse(Family == "", Order, 
                                         ifelse(Genus == "", Family, 
                                                ifelse(Species == "", 
                                                       Genus, paste(Genus, 
                                                                    Species, sep =" ")
                                                       )))))))
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

split_point <- 5000
# if <2500 ASVs, run in one set, else split and use a loop
# else run only once
# 2500 is safer with 8 Gb RAM but otherwise you can set this to 5000

split <- F
if(dim(seqtab.nochim)[2]>=split_point) split <-T


# set options for taxonomy assignment -------------------------------------

# check if this matches your dataset
pend <- paired_end # check if T or F
ovlp <- overlapping

if(!split){
  cat("assigning taxonomy...", "\n")
  
  taxtab_list <- Assign_taxonomy_ITS(seqtab = seqtab.nochim, 
                                 paired_end = pend, overlapping = ovlp)
  if(mean(is.na(taxtab_list$taxtab_slot[,2]))>0.2){
    RC<<-T
    taxtab_list <- Assign_taxonomy_ITS(seqtab = seqtab_nochim, 
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
    taxtab_list_temp <- Assign_taxonomy_ITS(seqtab = seqtab_nochim_temp, 
                                        paired_end = pend, overlapping = ovlp, rc = RC)
    if(i == 1 && mean(is.na(taxtab_list_temp$taxtab_slot[,2]))>0.2){
        RC<<-T
        taxtab_list_temp <- Assign_taxonomy_ITS(seqtab = seqtab_nochim_temp, 
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
beep(sound=sound_n)  
}

if(keep_time) toc()

if(play_audio) beep(sound=sound_n)

# fixing the output

taxtab <- taxtab_list$taxtab_slot
taxtab2 <- taxtab_list$taxtab2_slot

if(!split) {
  rm(taxtab_list, Assign_taxonomy_ITS)
} else {
  rm(taxtab_list, taxtab_list_temp, Assign_taxonomy_ITS) 
}
gc()

# fixing bad taxa ---------------------------------------------------------


# remove bad ides  ----------------------------------------

# which is the proportion of sequences with Kingdom only?
nASVs <- nrow(taxtab2)
nkingdom <- taxtab2 %>% dplyr::filter(Phylum == "") %>% nrow()
(f_nkingdom <-  nkingdom/nASVs)

# change as appropriate
filter_ASVs <- F
if(filter_ASVs){
  taxtab2_old <- taxtab2
  seqtab.nochim_old <- seqtab.nochim
  taxtab2 <- taxtab2 %>%
    dplyr::filter(Phylum != "")
  seqtab.nochim <- seqtab.nochim[,(colnames(seqtab.nochim) %in% taxtab2$ASV)]
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

dotree <- F # avoid if overlapping == F or if you have used the mixed option for merging
if(dim(seqtab.nochim)[2]>10000 | overlapping == F | merge_option == "mixed") {
  dotree <- F
}

if (dotree) {
  if(keep_time) tic("\nbuilding the phylogenetic tree")
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
  if(keep_time) toc()
}
if(play_audio) beep(sound = sound_n)

save.image(file = str_c(Study,"_small.Rdata"))

# Combine data into a phyloseq object ------------------------------------- 

metadata_path <- file.path("data", "metadata", "SraRunTable.txt") 
# alternatives when using data downloades from NCBI SRA are:
# "data/metadata/SraRunInfo.txt" 
# "data/metadata/SraRunTable.txt.csv"
# "data/metadata/SraRunTable.txt"
samdf <- read_tsv(metadata_path) 
if(ncol(samdf)==1) samdf <- read_csv(metadata_path)
# have a look at the sample data:
# if both bacteria and yeasts are detected, you should generate an accession list for filtering
# the following is an example of code for creating accession lists

run_acc_list_code <- F
if(run_acc_list_code){
  acc_list_bacteria <- samdf |>
    arrange(Run) |>
    dplyr::select(Run)
  write_tsv(acc_list_bacteria, file = "acc_list_bacteria.txt")
  acc_list_fungi <- samdf |>
    arrange(Run) |>
    dplyr::select(Run)
  write_tsv(acc_list_fungi, file = "acc_list_fungi.txt")
}

# use accession list to control which sequences will be processed
# only sequences in the accn_list file will be processed
use_accn_list <- F
accn_list_name <- "acc_list_fungi.txt" # need to adapt this to your own file
if(use_accn_list) {
  accn_list <- pull(read_tsv(accn_list_name, col_names = F),1)
  seq_to_process <- sample.names %in% accn_list
  fnFs <- fnFs[seq_to_process]
  if(paired_end) fnRs <- fnRs[seq_to_process]
  sample.names <- sample.names[seq_to_process]
}

if(use_accn_list) samdf <- dplyr::filter(samdf, Run %in% sample.names)

if(data_type == "sra"){
  if(all(rownames(seqtab.nochim) %in% samdf$Run)){
    cat("\nsamples in fastq files match samples in metadata\n")
  } else {
    cat("\nWARNING samples in fastq files DO NOT match samples in metadata\n")
  }
} else {
  if(all(rownames(seqtab.nochim) %in% samdf$Library_Name)){
    cat("\nsamples in fastq files match samples in metadata\n")
  } else {
    cat("\nWARNING samples in fastq files DO NOT match samples in metadata\n")
  }
}

# check if any cols in seqtab have 0 sums
# check if any cols in seqtab have 0 sums
seq_sums <- rowSums(seqtab.nochim)
runs_to_keep <- which(seq_sums>0)
runs <- rownames(seqtab.nochim)[runs_to_keep]
if(data_type == "sra") {
  if(length(runs_to_keep) < dim(seqtab.nochim)[1]) {
    seqtab.nochim <- seqtab.nochim[runs_to_keep,]
    samdf <- dplyr::filter(samdf, Run %in% runs)
  }
} else {
  if(length(runs_to_keep) < dim(seqtab.nochim)[1]) {
    seqtab.nochim <- seqtab.nochim[runs_to_keep,]
    samdf <- dplyr::filter(samdf, Library_Name %in% runs)
  }
}

# sanity check
# sanity check
if(data_type == "sra") {
  all(rownames(seqtab.nochim) %in% samdf$Run)
} else {
  all(rownames(seqtab.nochim) %in% samdf$Library_Name)
} # must be TRUE
# rownames(seqtab.nochim) <- sapply(strsplit(rownames(seqtab.nochim), "_", fixed = T), `[`, 1)
# all(rownames(seqtab.nochim) %in% samdf$Run) # TRUE
# add final number of sequences and number of issues from track2
# add number of sequences to track2
seqs <- as.data.frame(rowSums(seqtab.nochim))
seqs <- rownames_to_column(seqs, var = "Run")
colnames(seqs)[2] <-"seqs"
track2 <- left_join(track2, select(data = seqs, sample_name = Run, seqs))

if(data_type == "sra") {
  samdf <- left_join(samdf, select(track2, Run = sample_name, n_reads2 = nonchim, 
                                   n_issues))
} else {
  samdf <- left_join(samdf, select(track2, Library_Name = sample_name, n_reads2 = nonchim, 
                                   n_issues))
}
samdf <- as.data.frame(samdf)
rownames(samdf) <- samdf$Run

track2 <- mutate(track2, prop_tax_filt = data/nonchim)

rownames(samdf) <- samdf$Run

# a transposed sequence table
seqtab_t <- t(seqtab.nochim)
dim(seqtab_t)

# may have dropped a few runs because they lacked seqs
if(data_type == "sra"){
  samdf_2 <- samdf %>% dplyr::filter(Run %in% colnames(seqtab_t))
} else {
  samdf_2 <- samdf %>% dplyr::filter(Library_Name %in% colnames(seqtab_t))
}

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
                tax_database = "UNITE", Seq_accn = seq_accn,
                samples = n_samples, DOI = DOI, geoloc = loc_list, 
                primer_f, primer_r, overlapping, paired_end)
# saves study info
write_tsv(study, str_c(Study,"_study.txt"))


# format sample information for FMBN -------------------------------------------

# this needs to be adapted for each study
# check naming of the geoloc info

samples <- samples %>%
  mutate(description = str_c("Wine source tracking", Library.Name, Sample.Name, sep =", "))
samples <- samples %>%
  mutate(Sample_Name = Run) 

# information of geoloc (and names of the field) is very inconsistent:
# check the info in your sample metadata and adatp these commands
# use these if part or all of the geolocation information is missing
# samples$geo_loc_name_country <- "your country here"
# samples$geo_loc_name_country_continent <- "your continent here"
# samples$lat_lon <- NA_character_


# create label2 (to avoid numbers as first char.; s. can be removed later with
# tidyr::separate)
# need to be adjusted ad hoc especially for lat_lon
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
           geo_loc_country = geo_loc_name_country, 
           geo_loc_continent = geo_loc_name_country_continent, lat_lon, Sample_Name)
}

# save the sample information
write_tsv(samples, str_c(Study,"_samples.txt"))

# extract unique taxa -----------------------------------------------------

# a tibble with unique elements for taxonomy
# some changes are necessary for coherence with FMBN taxonomy
unique_tax <- taxtab2 %>% select(-ASV) %>% distinct()


# if the tax database is silva v138 class changes are not needed

unique_tax <- taxtab2 %>% select(-ASV) %>% distinct()

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


# create a phyloseq object FMBN style -------------------------------------

FMBN_physeq_OTU <- as.matrix(column_to_rownames(seqtab_ags, var = "s_label"))
FMBN_physeq_taxa <- as.matrix(column_to_rownames(unique_tax, var = "s_label"))
FMBN_physeq <- phyloseq(otu_table(FMBN_physeq_OTU, taxa_are_rows = T), 
                        sample_data(samples),
                        tax_table(FMBN_physeq_taxa))
saveRDS(FMBN_physeq, str_c("data/", Study, "_FMBN_ps.rds"))

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

#  save the log -----------------------------------------------------------
# need to put it at the end of the script because if the script is used
# on several days, yout ned to reinitialize the log
if(use_logr){
  log_path <- file.path(str_c(Study,"log", sep = "_"))
  log_open(log_path)
  gen_options <- list(
    audio = play_audio,
    sound = sound_n,
    timing = keep_time,
    verbose = verbose_output
  )
  log_print(gen_options, console = F)
  study_options <- list(
    study = Study,
    target_region <- str_c(target, region, sep = ", "),
    doi = DOI,
    Platform = platform
  )
  log_print(study_options, console = F)
  primers <- list(primerf = primer_f,
                  primerr = primer_r)
  log_print(primers, console = F)
  dotree
  log_close()
}

# Package citations -------------------------------------------------------

map(c(.cran_packages, .bioc_packages), citation)



# Credits and copyright ---------------------------------------------------

# Most of the script is taken from https://benjjneb.github.io/dada2/tutorial.html
# or https://benjjneb.github.io/dada2/bigdata.html
# with some changes and adaptations

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



