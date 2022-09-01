#########1#########2#########3#########4#########5#########6#########7#########8
# lactobacillaceae_ref_seqs v1.0 5/8//2022
#########1#########2#########3#########4#########5#########6#########7#########8

# description -------------------------------------------------------------

# NOTE FOR SELF can I find a way of using the tidysq package?
# NOTE FOR SELF can I use in some way the usedist package?
# NOTE FOR SELF Qiao et al., 2022 indicates 98.65% as a threshold for new species
# based on 16S RNA similarity, but only in terms of an isolate being a new
# species if similarity to any known species is <98.65
# Qiao, N., Wittouck, S., Mattarelli, P., Zheng, J., Lebeer, S., Felis, G.E., Gänzle, M.G., 2022. After the storm—Perspectives on the taxonomy of Lactobacillaceae. Jds Commun. https://doi.org/10.3168/jdsc.2021-0183


# This script is designed to:
# 1. assemble a database of 16S RNA gene reference sequences, including both 
#    sequences for type strains and reference sequences from the SILVA reference database
# 2. evaluate the ability to correctly assign type strain sequences (full 
#    length and V1-V3, V3-V4 and V4) to genus and species


# To run the script you need (in the project folder or working directory)
#   1. the metadata for type strains as a .tab delimited file; here 
#      Lactobacillaceae_metadata_w_outg.txt (it includes and outgroup column)
#   2. the sequences of the 16S RNA gene for the type strain, as a list lb_ts_seqs.rds
#   3. tab delimited file with data frames of reference sequences from SILVA
#      SILVA_genus_species.txt This is too large to upload on GitHub but can be 
#      easily created using the make_seq_df.R script and the SILVA v138.1 database
#      available on Zenodo (https://doi.org/10.5281/zenodo.4587954)
#   4. the reference databases for genus and species assignment formatted for dada2
#      (https://benjjneb.github.io/dada2/training.html) in a folder called tax_db 
#      placed one level up from the project folder or working directory (not included)



# Install/load packages ---------------------------------------------------

# may be buggy especially if there are newer versions of a package which
# need compilation

.cran_packages <- c("phylotools", "tidyverse", "kmer", "parallel", "caret", "broom",  
                    "beepr", "tictoc", "crayon", "data.table", "reshape2")
.bioc_packages <- c("BiocManager", "phyloseq", "DECIPHER", "phangorn", "dada2",
                    "treeio", "tidytree", "ggtree", "ggtreeExtra")

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
sapply(c(.cran_packages,.bioc_packages), require, character.only = TRUE)

# detect the number of cores
nc <- parallel::detectCores(logical = F) # to detect physical cores in MacOS
# set seed for reproducibility
set.seed(100L)
# options for keeping time and playing sound at the end of some functions
keep_time <- T
play_sound <- T


# functions ---------------------------------------------------------------


#  functions --------------------------------------------------------------

# cluster_seqs --------------------------------------------------------------

# this function is taken almost verbatim from the Bioconductor workflow for 
# microbiome data analysis (doi 10.12688/f1000research.8986.2)
# it takes an alignment as an input
# optimization depends on the boot_option (although anova on ML of trees
# created with either method shows that there is no difference, bootstrap
# values will always be 100% with stochastic rearrangement)

cluster_seqs <- function(alignment, boot_option = F){
  if(keep_time) tic("\nbuilding the phylogenetic tree")
  # transform in phydat object
  phang.align <- phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
  # create distance matrix
  cat(red("\ncreating distance matrix...","\n"))
  dm <- phangorn::dist.ml(phang.align)
  # perform Neighbor joining
  cat(red("\ncreating tree...","\n"))
  treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order
  # internal maximum likelihood for tree
  cat(red("\nestimating internal ML for tree...", "\n"))
  fit = phangorn::pml(treeNJ, data = phang.align)
  Sys.sleep(5)
  fitGTR <- update(fit, k = 4, inv = 0.2)
  Sys.sleep(5)
  # this is the step taking the longest time
  cat(red("\noptimization, please be patient...","\n"))
  if(boot_option){
    fitGTR <- optim.pml(
      fitGTR,
      model = "GTR",
      optNni = T
    )
  } else {
    fitGTR <- optim.pml(
      fitGTR,
      model = "GTR",
      optInv = TRUE,
      optGamma = TRUE,
      rearrangement = "stochastic",
      control = pml.control(trace = 0)
    )
  }
  cat(red("\nfitGTR object created\n"))
  # this is actually an object/list which contains the tree in the slot $tree
  # and the alignment in slot $data
  
  if(keep_time) toc()
  if(play_sound) beep(sound = 6)
  return(fitGTR)
}


# boot_tree ---------------------------------------------------------------

# a function for bootstrapping maximum likelihood phylogenetic trees

boot_tree <- function(fit, n_boot = 100, mcore = T){
  bs_pml <- bootstrap.pml(fit, bs = n_boot, multicore = mcore, mc.cores = nc) # 
  treeBS <- plotBS(fit$tree, bs_pml)
  return(treeBS)
}

# make_confusion_matrix ---------------------------------------------------

# create a function for extracting a confusion matrix in a tidy format
# no error trapping, you must load caret and broom before running
# takes a data frame containing the Genus_SILVA and Genus columns and uses the second as
# reference for the confusion matrix; Carnobacterium is filtered out

make_conf_matrices <- function(df){
  genus_pred_vs_ref <- df %>%
    dplyr::filter(Genus != "Carnobacterium") %>%
    dplyr::select(Genus_SILVA, Genus) %>%
    mutate(across(everything(), factor))
  # unify the levels
  unif_levels <- fct_unify(genus_pred_vs_ref)
  # estimate the confusion matrix
  conf_mat <- caret::confusionMatrix(data = unif_levels$Genus_SILVA,
                                     reference = unif_levels$Genus)
  # tidy with broom
  conf_mat_tidy <- tidy(conf_mat)
  return(conf_mat_tidy)
}

# kmer analysis -----------------------------------------------------------

# This function
# a. takes as an input a named character vector with sequences (argument seq_char)
# b. converts it to DNA.bin
# c. uses functions of package kmer with kmer_n length (default = 5) to
#    1. optionally calculate kmer distance matrix (calculated with mbed optional seed_vector if kmer_n>5)
#    2. return OTU membership vector at two similarity levels (sim_lvl_1 = 0.97, sim_lvl_2 = 0.99)
#    3. optionally perform MDS and return results and coordinates (which otherwise are null)
# results ar returned as a list

# I am not doing any effort for at error trapping (not even checking if required packages
# have been loaded: anyway kmer is needed and also ape, as a dependency, but it is loaded
# by phylotools), and crayon, for messages. It might also worth adding an option
# for the method in kmer::otu(); right now defaults to farthest, which I think is OK


kmer_analysis <- function(seq_char, kmer_n = 5, 
                          sim_lvl_1 = 0.97, sim_lvl_2 = 0.99,
                          calc_dist = T, perform_MDS = T, seed_vector = NULL){
  # create DNA_bin object
  cat(red("\ncreating DNAbin object...","\n"))
  seq_DNAbin <- as.DNAbin(DNAStringSet(seq_char))
  # creating the distance matrix
  
  if(calc_dist) {
    cat(red("\ncreating kmer distance matrix...","\n"))  
    if(kmer_n>5){
      kdist_mat <- mbed(seq_DNAbin, k = kmer_n, seeds = seed_vector)
    } else {
      kdist_mat <- kdistance(seq_DNAbin, k = kmer_n) 
    }
  } else {
    kdist_mat = NULL
  }
  
  cat(red("\ncalculating OTUs, similarity lvl 1", sim_lvl_1, "...", "\n"))
  otus_lvl1 <- otu(seq_DNAbin, k = kmer_n, threshold = sim_lvl_1, method = "farthest", nstart = 20)
  
  otus_lvl1_df <- tibble::enframe(otus_lvl1) %>%
    rename(seq_abbr = name, OTU = value) %>%
    mutate(seq_abbr = str_remove(seq_abbr, "\\*"))
  
  cat(red("\ncalculating OTUs, similarity lvl 2", sim_lvl_2, "...", "\n"))
  otus_lvl2 <- otu(seq_DNAbin, k = kmer_n, threshold = sim_lvl_2, method = "farthest", nstart = 20)
  
  otus_lvl2_df <- tibble::enframe(otus_lvl2) %>%
    rename(seq_abbr = name, OTU = value) %>%
    mutate(seq_abbr = str_remove(seq_abbr, "\\*"))
  
  if(calc_dist & perform_MDS) {
    cat(red("\nMDS...", "\n"))
    # this is a PCOA
    MDS_result <- cmdscale(kdist_mat, list. =T)
    MDS_coord <- as.data.frame(MDS_result$points) %>%
      rownames_to_column(var = "seq_abbr") %>%
      rename(dim1 = V1, dim2 = V2)
  } else {
    MDS_result <- MDS_coord <- NULL
  }
  
  
  kmer_list <- list(
    kmer_length = kmer_n,
    OTUsimlvl1 = sim_lvl_1,
    OTUsimlvl2 = sim_lvl_2,
    DNAbin = seq_DNAbin,
    kmer_dist_mat = kdist_mat,
    otus_1 = otus_lvl1,
    otus_2 = otus_lvl2,
    otus_1_df = otus_lvl1_df,
    otus_2_df = otus_lvl2_df,
    MDS_res = MDS_result,
    MDS_coord_df = MDS_coord
  )
  
  return(kmer_list)
}



# set a few other options -------------------------------------------------

taxo_group <- "Lactobacillaceae" # the taxonomic group you are working on


# setup section ends here -------------------------------------------------

# load 16S sequences and metadata for type strains ------------------------
# the type strain sequences, as a list

ts_seq_list <- readRDS(file = "lb_ts_seqs.rds")

ts_metadata <- read_tsv("Lactobacillaceae_metadata_w_outg.txt")

outgroups <- ts_metadata %>% 
  dplyr::select(Genus, group_16S, group_color,outgroup) %>%
  distinct(Genus, group_16S, group_color,outgroup)

genera <- ts_metadata %>% 
  dplyr::filter(Genus != "Carnobacterium") %>%
  select(Genus) %>% 
  distinct() %>% pull()

# load and filter Silva 138.1 reference sequences ---------------------------

# the tab-delimited file was created using the script make_seqs_df

# only need the genus_species file, sequences are identical
SILVA_genus_species <- read_tsv("SILVA_genus_species.txt")
SILVA_genus_species <- SILVA_genus_species %>%
  dplyr::filter(Genus %in% genera)


# re-identification of type strain sequences with SILVA ----------------------------

# carry out re-identification of ts sequences using SILVA v138.1 for:
# a. full sequence
# b. V1-V3
# c. V3-V4
# d. V4
# all sequences are in the right orientation and none contains long stretches of N

# pointing to the taxonomic reference databases
taxdb_dir <- file.path("..","tax_db") # change this if the tax databases are elsewhere
ref_fasta <- file.path(taxdb_dir, "silva_nr99_v138_1_train_set.fa")
sp_ass <- file.path(taxdb_dir, "silva_species_assignment_v138_1.fa")

# create a list to accumulate the results
# careful here, the loop depends on the structure of this list
ts_ides <- vector(mode = "list", length = length(ts_seq_list))
# careful here, I am not attempting error trapping

if(keep_time) tic("\ntaxonomy assignment with SILVA")

for(i in seq_along(ts_seq_list)){
  # create a pretend seq table and remove duplicates
  seqs <- ts_seq_list[[i]][[1]]
  dupl_seqs <- seqs[which(duplicated(seqs))]
  seqs <- seqs[-which(duplicated(seqs))]
  pretend_seqtab <- as_tibble(seqs) %>% 
    rename(seq.text = value) %>%
    mutate(seq.name = names(seqs),
           dummy1 = 1,
           dummy2 = 1
    ) %>%
    column_to_rownames(var = "seq.text")
  cat(red("\nassigning taxonomy, genus, iteration", i, names(ts_seq_list)[i], "\n"))
  seq_matrix <- pretend_seqtab %>%
    dplyr::select(-seq.name)
  seq_matrix <- t(as.matrix(seq_matrix))
  seq_matrix_no_ambiguities <- seq_matrix
  colnames(seq_matrix_no_ambiguities) <- str_remove_all(colnames(seq_matrix_no_ambiguities),"[RYWSMKHBVDNSU]")
  
  taxtab <- assignTaxonomy(seq_matrix_no_ambiguities, 
                           refFasta = ref_fasta, multithread = TRUE)
  
  cat(red("\nassigning taxonomy, species, iteration", i, names(ts_seq_list)[i], "\n"))
  taxtab <- addSpecies(taxtab, refFasta = sp_ass)
  taxtab_df <- as.data.frame(taxtab) %>%
    rownames_to_column(var = "mod_seq") %>%
    bind_cols(rownames_to_column(pretend_seqtab, var = "seq.text")) %>%
    dplyr::rename(Genus_SILVA = Genus, Species_SILVA = Species) %>%
    dplyr::select(-dummy1, -dummy2)
  # add back the duplicates
  dupli_df <- data.frame(seq.text = dupl_seqs, dupli_strain = names(dupl_seqs)) %>%
    left_join(taxtab_df) %>% mutate(seq.name = dupli_strain) %>%
    select(mod_seq, Kingdom:Species_SILVA, seq.text, seq.name)
  # populate the list
  ts_ides[[i]] <- bind_rows(taxtab_df, dupli_df) %>% arrange(seq.name)
  names(ts_ides)[i] <- names(ts_seq_list)[i]
}
if(keep_time) toc()
if(play_sound) beep(sound = 6)
rm(seqs, dupl_seqs, pretend_seqtab, seq_matrix, seq_matrix_no_ambiguities, taxtab, taxtab_df, dupli_df)

# adding metadata
ts_ides_ann <- map(ts_ides, ~ left_join(.x, dplyr::select(ts_metadata, seq.name = seq_abbr, Genus:Species)))
# View the data frames
# View(ts_ides_ann[[1]])
# View(ts_ides_ann[[2]])
# View(ts_ides_ann[[3]])
# View(ts_ides_ann[[4]])


# create confusion matrices for taxo assignment ---------------------------

# use map for the confusion matrices
conf_matrices <- map_dfr(ts_ides_ann, make_conf_matrices, .id = "region")

conf_matrices %>% 
  dplyr::filter(!is.na(class) & (term == "sensitivity" | term == "specificity")) %>% 
  ggplot(mapping = aes(x = class, y = estimate, shape = region)) +
  facet_wrap(~term, scales = "free_y") +
  geom_point() +
  labs(x = "Genus") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

# ggsave(file = str_c(taxo_group,"_confusionmatrices.tiff"), dpi = 300)

# prepare reference sequences ---------------------------------------------------------

SILVA_genus_species <- SILVA_genus_species %>% 
  mutate(rs_id = str_c("rs",seq(1:nrow((SILVA_genus_species)))),
         sp_2 = if_else(is.na(Species), "",Species)) %>%
  mutate(seq_abbr = str_c(
    str_trunc(Genus, 7, side = "right", ellipsis = "_"),
    str_trunc(sp_2, 4, side = "right", ellipsis = "_"),
    rs_id
  )) %>%
  mutate(seq.name = seq_abbr) %>%
  select(-sp_2)

# saveRDS(SILVA_genus_species, file = "Lb_SILVA_genus_species_mdata.RDS")

# create vectors of sequences

ref_seq_fl <- SILVA_genus_species$seq.text
names(ref_seq_fl) <- SILVA_genus_species$seq.name

# extract regions
p28f <- "TTTGATCNTGGCTC"
p519R <- as.character(reverseComplement(DNAStringSet("GTNTTACNGCGGCKGCTG")))

p341f <- "CCTACGGGNGGCWGCAG"
p515F <- "GTGCCAGCMGCCGCGGTAA"
p805R <- as.character(reverseComplement(DNAStringSet("GACTACHVGGGTATCTAATCC")))

# NOTE FOR SELF a tidy output is obtained with tidysq::find_motif()
primer_matches_28f <- vmatchPattern(p28f, DNAStringSet(ref_seq_fl),
                                    max.mismatch = 2, min.mismatch = 0, fixed = F)
primer_matches_341f <- vmatchPattern(p341f, DNAStringSet(ref_seq_fl),
                                     max.mismatch = 2, min.mismatch = 0, fixed = F)
primer_matches_515F <- vmatchPattern(p515F, DNAStringSet(ref_seq_fl),
                                     max.mismatch = 2, min.mismatch = 0, fixed = F)
primer_matches_519R <- vmatchPattern(p519R, DNAStringSet(ref_seq_fl),
                                     max.mismatch = 2, min.mismatch = 0, fixed = F)
primer_matches_805R <- vmatchPattern(p805R, DNAStringSet(ref_seq_fl),
                                     max.mismatch = 2, min.mismatch = 0, fixed = F)

v1v3_end <- primer_matches_519R@ends
# some elements are of length >1
v1v3_end_2 <- vector(mode = "integer", length = length(v1v3_end))
for(i in seq_along(v1v3_end)){
  if(is.null(v1v3_end[[i]])) { 
    v1v3_end_2[[i]] <- NA_integer_
  } else {
    v1v3_end_2[[i]] <- v1v3_end[[i]][1]
  }
}

# v1v3_end[sapply(v1v3_end, is.null)]<-NA_integer_
v1v3_start_end <- data.frame(
  string = ref_seq_fl,
  start = rep(1, length(ref_seq_fl)),
  end = v1v3_end_2)


# any NA?
which(is.na(v1v3_start_end$end))
# ad hoc change
mean_end_V3 <- ceiling(mean(v1v3_start_end$end, na.rm = T))
v1v3_start_end$end[which(is.na(v1v3_start_end$end))] <- mean_end_V3
which(is.na(v1v3_start_end$end))

# now extract 
ref_seq_V1V3 <- pmap_chr(v1v3_start_end, str_sub)
names(ref_seq_V1V3) <- names(ref_seq_fl)

# same V3V4
v3v4_start <- primer_matches_341f@ends

v3v4_start_2 <- vector(mode = "integer", length = length(v3v4_start))
for(i in seq_along(v3v4_start)){
  if(is.null(v3v4_start[[i]])) { 
    v3v4_start[[i]] <- NA_integer_
  } else {
    v3v4_start_2[[i]] <- v3v4_start[[i]][1]
  }
}


v3v4_end <- primer_matches_805R@ends

# some elements are of length >1
v3v4_end_2 <- vector(mode = "integer", length = length(v3v4_end))
for(i in seq_along(v3v4_end)){
  if(is.null(v3v4_end[[i]])) { 
    v3v4_end_2[[i]] <- NA_integer_
  } else {
    v3v4_end_2[[i]] <- v3v4_end[[i]][1]
  }
}

v3v4_start_end <- data.frame(
  string = ref_seq_fl,
  start = v3v4_start_2,
  end = v3v4_end_2)


# patching up
mean_start_v3v4 <- ceiling(mean(v3v4_start_end$start, na.rm = T))
mean_end_v3v4 <- ceiling(mean(v3v4_start_end$end, na.rm = T))
v3v4_start_end$start[which(is.na(v3v4_start_end$start))] <- mean_start_v3v4
v3v4_start_end$end[which(is.na(v3v4_start_end$end))] <- mean_end_v3v4
which(is.na(v3v4_start_end$start))
which(is.na(v3v4_start_end$end))

# extracting sequences and aligning
ref_seq_V3V4 <- pmap_chr(v3v4_start_end, str_sub)
names(ref_seq_V3V4) <- names(ref_seq_fl)


# same for V4

v4_start <- primer_matches_515F@ends

v4_start_2 <- vector(mode = "integer", length = length(v4_start))
for(i in seq_along(v4_start)){
  if(is.null(v4_start[[i]])) { 
    v4_start[[i]] <- NA_integer_
  } else {
    v4_start_2[[i]] <- v4_start[[i]][1]
  }
}



v4_start_end <- data.frame(
  string = ref_seq_fl,
  start = v4_start_2,
  end = v3v4_end_2)


# patching up
mean_start_v4 <- ceiling(mean(v4_start_end$start, na.rm = T))
mean_end_v4 <- ceiling(mean(v4_start_end$end, na.rm = T))
v4_start_end$start[which(is.na(v4_start_end$start))] <- mean_start_v4
v4_start_end$end[which(is.na(v4_start_end$end))] <- mean_end_v4
which(is.na(v4_start_end$start))
which(is.na(v4_start_end$end))

# extracting sequences and aligning
ref_seq_V4 <- pmap_chr(v4_start_end, str_sub)
names(ref_seq_V4) <- names(ref_seq_fl)

# concatenating the sequence vectors

ref_seq_list <- list(
  full_length = c(ts_seq_list$full_sequence$seqs, ref_seq_fl),
  V1V3 = c(ts_seq_list$V1V3$seqs, ref_seq_V1V3),
  V3V4 = c(ts_seq_list$V3V4$seqs, ref_seq_V3V4),
  V4 = c(ts_seq_list$V1V3$seqs, ref_seq_V4)
)
save(ref_seq_list, file = str_c(taxo_group, "_ref_seqs_list.Rdata", sep = ""))

# perform kmer analysis ---------------------------------------------------

# perform kmer analysis on type strains and reference sequences 
# I am not bothering to create a loop or using a functional
# use 0.97 and 0.9865

# on the full length sequences
# the first 374 sequences are type strains
if (keep_time) tic("kmer analysis, full sequence")
kmer_fl_ts_refseq <- kmer_analysis(ref_seq_list$full_length, kmer_n = 5, 
                                   sim_lvl_2 = 0.9865,
                                   calc_dist = T,
                                   perform_MDS = F)
if (keep_time) toc()


# V1-V3
hist(nchar(ref_seq_list$V1V3))
summary(nchar(ref_seq_list$V1V3))
if (keep_time) tic("kmer analysis, V1V3")
kmer_V1V3_ts_refseq <- kmer_analysis(ref_seq_list$V1V3[nchar(ref_seq_list$V1V3)>400 & nchar(ref_seq_list$V1V3)<560], kmer_n = 5, 
                                     sim_lvl_2 = 0.9865,
                                     calc_dist = T,
                                     perform_MDS = F)
if (keep_time) toc()


# V3-V4
hist(nchar(ref_seq_list$V3V4))
summary(nchar(ref_seq_list$V3V4))
if (keep_time) tic("kmer analysis, V3V4")
kmer_V3V4_ts_refseq <- kmer_analysis(ref_seq_list$V3V4[nchar(ref_seq_list$V3V4)>300 & nchar(ref_seq_list$V3V4)<460], kmer_n = 5, 
                                     sim_lvl_2 = 0.9865,
                                     calc_dist = T,
                                     perform_MDS = F)
if (keep_time) toc()

# V4
hist(nchar(ref_seq_list$V4))
summary(nchar(ref_seq_list$V4))
if (keep_time) tic("kmer analysis, V4")
kmer_V4_ts_refseq <- kmer_analysis(ref_seq_list$V4[nchar(ref_seq_list$V4)>200 & nchar(ref_seq_list$V4)>400], kmer_n = 5, 
                                   calc_dist = T,
                                   sim_lvl_2 = 0.9865,
                                   perform_MDS = F)
if (keep_time) toc()

n_distinct(kmer_fl_ts_refseq$otus_2)

# too many OTUs, I need to find a distance based criterion

# library(reshape2)
dist_df_fl <- reshape2::melt(as.matrix(kmer_fl_ts_refseq$kmer_dist_mat))
# remove 0 (diagonal)
# only keep type strains in the second column and reference in the first
ts_abbr <- ts_metadata %>% pull(seq_abbr)
ref_abbr <- SILVA_genus_species %>% pull(seq_abbr)
dist_df_fl <- dist_df_fl %>% dplyr::filter(value != 0) %>%
  dplyr::filter(Var1 %in% ref_abbr) %>%
  dplyr::filter(Var2 %in% ts_abbr) %>%
  arrange(Var1, value)

qplot(value, data = dist_df_fl, geom = "density")

# keep minimum distance
min_dist_fl <- dist_df_fl %>% group_by(Var1) %>% 
  summarise(value = min(value)) %>% left_join(dist_df_fl)
qplot(value, data = min_dist_fl, geom = "density")
# remove duplicates
min_dist_fl <- min_dist_fl %>% dplyr::filter(!duplicated(min_dist_fl$Var1))

# merge taxonomic information

min_dist_fl <- left_join(min_dist_fl, select(ts_metadata, seq_abbr, group_16S, Genus, Species),
                         by = c("Var2" = "seq_abbr")) 
# add further taxonomic information
min_dist_fl <- min_dist_fl %>%
  mutate(Kingdom = "Bacteria",
         Phylum = "Firmicutes",
         Class = "Bacilli",
         Order = "Lactobacillales", 
         Family = "Lactobacillaceae")

# add sequences 
min_dist_fl <- left_join(min_dist_fl, select(SILVA_genus_species, Var1 = seq_abbr, seq.text))


# create and save custom species reference ----------------------------------------

custom_species_ref_df <- select(min_dist_fl, Var1, Genus, Species, seq.text) %>%
  tidyr::unite("seq.name", Var1, Genus, Species, sep = " ")
# let's add the type strains
type_strains_seq_df <- select(ts_ides_ann$full_sequence, seq.text = mod_seq, seq_abbr = seq.name, Genus, Species)
# add strain info and create names
type_strains_seq_df <- type_strains_seq_df %>% 
  left_join(select(ts_metadata, seq_abbr, Accn_n)) %>%
  tidyr::unite("seq.name", Accn_n, Genus, Species, sep = " ") %>%
  select(seq.name, seq.text)
custom_species_ref_df <-bind_rows(custom_species_ref_df, type_strains_seq_df)
phylotools::dat2fasta(custom_species_ref_df, outfile = "Lactobacillaceae_custom_sp_ass.fa")

# try species assignment with custom reference database -------------------

sp_ass <- "Lactobacillaceae_custom_sp_ass.fa"
ts_ides_2 <- vector(mode = "list", length = length(ts_seq_list))
# careful here, I am not attempting error trapping
if(keep_time) tic("\ntaxonomy assignment with SILVA")
for(i in seq_along(ts_seq_list)){
  seqs <- ts_seq_list[[i]][[1]]
  dupl_seqs <- seqs[which(duplicated(seqs))]
  seqs <- seqs[-which(duplicated(seqs))]
  seqtab <- as_tibble(seqs) %>% 
    mutate(seq.name = names(seqs)) %>%
    rename(seq.text = value)
  seqs<-str_remove_all(seqs,"[RYWSMKHBVDNSU]")
  cat(red("\nassigning taxonomy, genus, iteration", i, names(ts_seq_list)[i], "\n"))
  taxtab <- assignTaxonomy(seqs, 
                           refFasta = ref_fasta, multithread = TRUE)
  cat(red("\nassigning taxonomy, species, iteration", i, names(ts_seq_list)[i], "\n"))
  taxtab <- addSpecies(taxtab, refFasta = sp_ass, allowMultiple = 3)
  taxtab_df <- as.data.frame(taxtab) %>%
    rownames_to_column(var = "mod_seq") %>%
    bind_cols(seqtab) %>%
    dplyr::rename(Genus_SILVA = Genus, Species_SILVA = Species) 
  # add back the duplicates
  dupli_df <- data.frame(seq.text = dupl_seqs, dupli_strain = names(dupl_seqs)) %>%
    left_join(taxtab_df) %>% mutate(seq.name = dupli_strain) %>%
    select(mod_seq, Kingdom:Species_SILVA, seq.text, seq.name)
  # populate the list
  ts_ides_2[[i]] <- bind_rows(taxtab_df, dupli_df) %>% arrange(seq.name)
  names(ts_ides_2)[i] <- names(ts_seq_list)[i]
}
if(keep_time) toc()
if(play_sound) beep(sound = 6)
rm(seqs, dupl_seqs, taxtab, seqtab, taxtab_df, dupli_df)

# adding metadata
ts_ides_2_ann <- map(ts_ides_2, ~ left_join(.x, dplyr::select(ts_metadata, seq.name = seq_abbr, Genus:Species)))
View(ts_ides_2_ann[[1]])
View(ts_ides_2_ann[[2]])
View(ts_ides_2_ann[[3]])
View(ts_ides_2_ann[[4]])

# assemble in data frames (I am not bothering to build a function and use map)

sp_matches_ts_fl <- ts_ides_2_ann[[1]] %>% 
  mutate(species_match = case_when(
    Species_SILVA == Species ~ "match",
    str_detect(Species_SILVA, "/") & str_detect(Species_SILVA, Species) ~ "pmatch",
    TRUE ~ "nomatch"
  )) %>% 
  dplyr::filter(Genus != "Carnobacterium") %>%
  group_by(Genus) %>% 
  count(species_match) %>%
  pivot_wider(names_from = species_match, values_from = n, values_fill = 0) %>%
  mutate(prop_matches = (match/(match+pmatch+nomatch)))


sp_matches_ts_V1V3 <- ts_ides_2_ann[[2]] %>% 
  mutate(species_match = case_when(
    Species_SILVA == Species ~ "match",
    str_detect(Species_SILVA, "/") & str_detect(Species_SILVA, Species) ~ "pmatch",
    TRUE ~ "nomatch"
  )) %>% 
  dplyr::filter(Genus != "Carnobacterium") %>%
  group_by(Genus) %>% 
  count(species_match) %>%
  pivot_wider(names_from = species_match, values_from = n, values_fill = 0) %>%
  mutate(prop_matches = (match/(match+pmatch+nomatch)))

sp_matches_ts_V3V4 <- ts_ides_2_ann[[3]] %>% 
  mutate(species_match = case_when(
    Species_SILVA == Species ~ "match",
    str_detect(Species_SILVA, "/") & str_detect(Species_SILVA, Species) ~ "pmatch",
    TRUE ~ "nomatch"
  )) %>% 
  dplyr::filter(Genus != "Carnobacterium") %>%
  group_by(Genus) %>% 
  count(species_match) %>%
  pivot_wider(names_from = species_match, values_from = n, values_fill = 0) %>%
  mutate(prop_matches = (match/(match+pmatch+nomatch)))


sp_matches_ts_V4 <- ts_ides_2_ann[[4]] %>% 
  mutate(species_match = case_when(
    Species_SILVA == Species ~ "match",
    str_detect(Species_SILVA, "/") & str_detect(Species_SILVA, Species) ~ "pmatch",
    TRUE ~ "nomatch"
  )) %>% 
  dplyr::filter(Genus != "Carnobacterium") %>%
  group_by(Genus) %>% 
  count(species_match) %>%
  pivot_wider(names_from = species_match, values_from = n, values_fill = 0) %>%
  mutate(prop_matches = (match/(match+pmatch+nomatch)))

summary(sp_matches_ts_fl)

fl <- sp_matches_ts_fl %>% 
  ungroup() %>% 
  summarize(no_match = sum(nomatch),
            match = sum(match),
            pmatch = sum(pmatch))
V1V3 <- sp_matches_ts_V1V3 %>% 
  ungroup() %>% 
  summarize(no_match = sum(nomatch),
            match = sum(match),
            pmatch = sum(pmatch))
V3V4 <- sp_matches_ts_V3V4 %>% 
  ungroup() %>% 
  summarize(no_match = sum(nomatch),
            match = sum(match),
            pmatch = sum(pmatch))
V4 <- sp_matches_ts_V4 %>% 
  ungroup() %>% 
  summarize(no_match = sum(nomatch),
            match = sum(match),
            pmatch = sum(pmatch))

sp_matches <- bind_rows(fl, V1V3, V3V4, V4)
sp_matches$region <- c("full length", "V1V3", "V3V4", "V4")

# write_tsv(sp_matches, file = str_c(taxo_group, "_typestrain_spmatches.txt"))

# create bar plots

match_df <- bind_rows(
  mutate(sp_matches_ts_fl, region = "full length"),
  mutate(sp_matches_ts_V1V3, region = "V1-V3"),
  mutate(sp_matches_ts_V3V4, region = "V3-V4"),
  mutate(sp_matches_ts_V4, region = "V4"),
) %>% select(-prop_matches) %>%
  pivot_longer(cols = nomatch:pmatch,
               names_to = "match_type",
               values_to = "weight") %>%
  mutate(match_type = factor(match_type, levels = c("nomatch", "pmatch", "match"))) %>%
  ungroup()

sp_match_bar_plot <- ggplot(match_df, aes(x = Genus, fill = match_type, y = weight)) + 
  facet_wrap(~region)+
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("nomatch" = "red", "pmatch" = "cyan", "match" = "blue")) +
  labs(fill = "species assignm.", y = "proportion") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
sp_match_bar_plot


#ggsave(sp_match_bar_plot, filename = str_c(taxo_group,"spmatchesbarplot.tiff", sep = "_"), dpi = 300)

match_sp_ide_box_plot <- bind_rows(
  mutate(sp_matches_ts_fl, region = "full length"),
  mutate(sp_matches_ts_V1V3, region = "V1-V3"),
  mutate(sp_matches_ts_V3V4, region = "V3-V4"),
  mutate(sp_matches_ts_V4, region = "V4") 
) %>% 
  ggplot(mapping = aes(x= region, y = prop_matches)) +
  geom_boxplot(notch = T) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(y = "proportion of matching species identification")
match_sp_ide_box_plot
ggsave(match_sp_ide_box_plot, filename = str_c(taxo_group,"spmatchesboxplot.tiff", sep = "_"), dpi = 300)
  
# summary info on proportion of matches
bind_rows(
  mutate(sp_matches_ts_fl, region = "full length"),
  mutate(sp_matches_ts_V1V3, region = "V1-V3"),
  mutate(sp_matches_ts_V3V4, region = "V3-V4"),
  mutate(sp_matches_ts_V4, region = "V4"),
) %>% 
  ungroup() %>%
  group_by(region) %>%
  summarise (median_prop_match = median(prop_matches))



