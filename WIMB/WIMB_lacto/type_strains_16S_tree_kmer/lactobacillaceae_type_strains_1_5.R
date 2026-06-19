#########1#########2#########3#########4#########5#########6#########7#########8
# lactobacillaceae_type_strains v1_5 16/12/2022
#########1#########2#########3#########4#########5#########6#########7#########8
# This script is designed load a data frame containing 16S RNA gene sequences for
# type strains belonging to family, with further metadata
# and to create dendrograms and other representations of distance 
# among strains and species.
# Sequences have been manually downloaded from a variety of sources 
# (mainly LPSN, GenBank, RefSeq). 
# In addition, this script will create a dataframe with reference sequences
# obtained from the SILVA v138.1 reference database

# To run the script you need (in the project folder or working directory)
#   the data frame with metadata (ts_metadata_small.txt)
#   a fasta file with the sequences for the 16S RNA gene of the type strains 
#     (lactob_type_16S.fasta) (for convenience the aligned sequences and
#     a list containing the sequences for the full length gene and hypervariable
#     regions estracted from them are also provided in this example)
#   a R list with the sequences as named vectors (optionally)


# BEWARE: this is designed to run on UNIX like systems; processors should be set to NULL
# in Windows

# Install/load packages ---------------------------------------------------

# may be buggy especially if there are newer versions of a package which
# need compilation

.cran_packages <- c("phylotools", "tidyverse", "randomcoloR", "readxl", 
                    "RColorBrewer", "kmer", "vegan", "insect",
                    "parallel", "beepr", "tictoc", "crayon")
.bioc_packages <- c("BiocManager","dada2", "phyloseq", "DECIPHER", "phangorn", 
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
set.seed(100L) # if you want to perform bootstrap maybe use set.seed(NULL)
# options for keeping time and playing sound at the end of some functions
keep_time <- T
play_sound <- T
# option for loading fasta files; set to T if you want to rebuild the
# data frame for sequences, F otherwise; useful if new sequences are added
run_extra <- F


# set a few other options -------------------------------------------------

taxo_group <- "Lactobacillaceae" # the taxonomic group you are working on


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
    set.seed(NULL)
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

boot_tree <- function(fit, n_boot = 100, mcore = T, psupport = 80){
  bs_pml <- bootstrap.pml(fit, bs = n_boot, multicore = mcore, mc.cores = nc, 
                          optiNni = T, control = pml.control(trace = 0))  
  treeBS <- plotBS(tree = fit$tree, BStrees = bs_pml, p = psupport)
  return(treeBS)
}


# load fasta for type strains in a data frame ------------------------------

## Run only if you want to add new sequences
# a fast way to load multiple fasta files and attempt to parse the sequence information
# this works, for example, with fasta files downloaded from LPSN
if(run_extra){
  # put the path with your folder containing your fasta files in this command
  lactob_fasta_files <- list.files(file.path(".","fasta_Lactobacillaceae"))
  lactob_fasta_files_path <- file.path(".", "fasta_Lactobacillaceae", lactob_fasta_files)
  # read fasta files in a data frame
  lactob_ts_df <- map_dfr(.x = lactob_fasta_files_path, .f = phylotools::read.fasta, clean_name = F)
  
  lactob_ts_df <- lactob_ts_df %>% mutate(LPSN_seq_length = str_length(seq.text))
  
  strain_matrix <- str_split_fixed(lactob_ts_df$seq.name, "__", 3)
  colnames(strain_matrix) <- c("s_label", "Strain", "Accn_n")
  # now I need to parse species name and strain name
  species_strain_df <- as_tibble(strain_matrix) %>% 
    separate(s_label, into = c("Genus", "Species", "nothing", "subspecies"), sep = "_", remove = F) %>%
    select(-nothing)
  # may need further manual fixes for subspecies
  # you need to rebuild fasta file and metadata suing these
}

## end not run


# load data ---------------------------------------------------------------

# load the sequences from fasta
seqsdf <- phylotools::read.fasta(file = "lactob_type_16S.fasta")
seq_metadata <- read_tsv(file = "ts_metadata_small.txt")

# if you want to load the .rds object with sequences (which is faster)
# ts_seqs <- readRDS("lb_ts_seqs.rds")

# tidysq requires more time than phylotools but objects are smaller
# tidysq::read_fasta(file_name = "lactob_type_16S.fasta")
# naming of columns is also different

# create a character vector to be used in alignment

lb_type_16S_seqs <- seqsdf$seq.text
names(lb_type_16S_seqs) <- seqsdf$seq.name

#  create alignment -------------------------------------------------------

lactob_type_16S_aligned <-  DECIPHER::AlignSeqs(DNAStringSet(lb_type_16S_seqs), anchor = NA, processors = nc)

BrowseSeqs(lactob_type_16S_aligned)

writeXStringSet(lactob_type_16S_aligned, filepath = "lactob_type_16S_aligned.fasta", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")
# save the alignment
phangorn::write.phyDat(lactob_type_16S_aligned, file = "lactob_type_16S_aligned.fasta", format = "fasta")
# sequences vary in length
hist(nchar(lb_type_16S_seqs))

#  extracting sequences for regions and aligning --------------------------

# extracting regions ------------------------------------------------------

# frequently used primers are
# V1-V3 Gray28f	TTTGATCNTGGCTC Gray519r	GTNTTACNGCGGCKGCTG
# V3-V4 Bakt_341F	CCTACGGGNGGCWGCAG Bakt_805R	GACTACHVGGGTATCTAATCC
# V4 515F	GTGCCAGCMGCCGCGGTAA Bakt_805R	GACTACHVGGGTATCTAATCC

# this is to detect the start positions

p28f <- "TTTGATCNTGGCTC"
p519R <- as.character(reverseComplement(DNAStringSet("GTNTTACNGCGGCKGCTG")))

p341f <- "CCTACGGGNGGCWGCAG"
p515F <- "GTGCCAGCMGCCGCGGTAA"
p805R <- as.character(reverseComplement(DNAStringSet("GACTACHVGGGTATCTAATCC")))

# NOTE FOR SELF a tidy output is obtained with tidysq::find_motif()
primer_matches_28f <- vmatchPattern(p28f, DNAStringSet(lb_type_16S_seqs),
                                    max.mismatch = 2, min.mismatch = 0, fixed = F)
primer_matches_341f <- vmatchPattern(p341f, DNAStringSet(lb_type_16S_seqs),
                                     max.mismatch = 2, min.mismatch = 0, fixed = F)
primer_matches_515F <- vmatchPattern(p515F, DNAStringSet(lb_type_16S_seqs),
                                     max.mismatch = 2, min.mismatch = 0, fixed = F)
primer_matches_519R <- vmatchPattern(p519R, DNAStringSet(lb_type_16S_seqs),
                                     max.mismatch = 2, min.mismatch = 0, fixed = F)
primer_matches_805R <- vmatchPattern(p805R, DNAStringSet(lb_type_16S_seqs),
                                     max.mismatch = 2, min.mismatch = 0, fixed = F)

v1v3_end <- primer_matches_519R@ends
v1v3_end[sapply(v1v3_end, is.null)]<-NA_integer_
v1v3_start_end <- data.frame(
  string = lb_type_16S_seqs,
  start = rep(1, length(lb_type_16S_seqs)),
  end = base::unlist(v1v3_end, recursive = F)
)

# any NA?
which(is.na(v1v3_start_end$end))
# ad hoc change
mean_end_V3 <- ceiling(mean(v1v3_start_end$end, na.rm = T))
v1v3_start_end$end[which(is.na(v1v3_start_end$end))] <- mean_end_V3
which(is.na(v1v3_start_end$end))

# now extract 
lb_type_16S_seqs_V1V3 <- pmap_chr(v1v3_start_end, str_sub)
names(lb_type_16S_seqs_V1V3) <- names(lb_type_16S_seqs)
histogram(nchar(lb_type_16S_seqs_V1V3))
# one sequence is very long
which(nchar(lb_type_16S_seqs_V1V3) == max(nchar(lb_type_16S_seqs_V1V3)))


# now V3-V4
v3v4_start <- primer_matches_341f@ends
v3v4_start[sapply(v3v4_start, is.null)]<-NA_integer_
v3v4_end <- primer_matches_805R@ends
v3v4_end[sapply(v3v4_end, is.null)]<-NA_integer_
v3v4_start_end <- data.frame(
  string = lb_type_16S_seqs,
  start = base::unlist(v3v4_start, recursive = F),
  end = base::unlist(v3v4_end, recursive = F)
)

# patching up
mean_start_v3v4 <- ceiling(mean(v3v4_start_end$start, na.rm = T))
mean_end_v3v4 <- ceiling(mean(v3v4_start_end$end, na.rm = T))
v3v4_start_end$start[which(is.na(v3v4_start_end$start))] <- mean_start_v3v4
v3v4_start_end$end[which(is.na(v3v4_start_end$end))] <- mean_end_v3v4
which(is.na(v3v4_start_end$start))
which(is.na(v3v4_start_end$end))

# extracting sequences
lb_type_16S_seqs_V3V4 <- pmap_chr(v3v4_start_end, str_sub)
names(lb_type_16S_seqs_V3V4) <- names(lb_type_16S_seqs)
histogram(nchar(lb_type_16S_seqs_V3V4))


# finally, V4
v4_start <- primer_matches_515F@ends
v4_start[sapply(v4_start, is.null)]<-NA_integer_
v4_end <- primer_matches_805R@ends
v4_end[sapply(v4_end, is.null)]<-NA_integer_
v4_start_end <- data.frame(
  string = lb_type_16S_seqs,
  start = base::unlist(v4_start, recursive = F),
  end = base::unlist(v4_end, recursive = F)
)

# patching up
mean_start_v4 <- ceiling(mean(v4_start_end$start, na.rm = T))
mean_end_v4 <- ceiling(mean(v4_start_end$end, na.rm = T))
v4_start_end$start[which(is.na(v4_start_end$start))] <- mean_start_v4
v4_start_end$end[which(is.na(v4_start_end$end))] <- mean_end_v4
which(is.na(v4_start_end$start))
which(is.na(v4_start_end$end))

# extracting sequences
lb_type_16S_seqs_V4 <- pmap_chr(v4_start_end, str_sub)
names(lb_type_16S_seqs_V4) <- names(lb_type_16S_seqs)
histogram(nchar(lb_type_16S_seqs_V4))

# saving character vectors for sequences as a list
lb_ts_seqs <- list(
  FL = lb_type_16S_seqs,
  V1V3 = lb_type_16S_seqs_V1V3,
  V3V4 = lb_type_16S_seqs_V3V4,
  V4 = lb_type_16S_seqs_V4
)
saveRDS(lb_ts_seqs, file = "lb_ts_seqs.rds")

# align hypervariable regions ---------------------------------------------

# V1-V3
lactob_type_16S_V1V3_aligned <- DECIPHER::AlignSeqs(DNAStringSet(lb_type_16S_seqs_V1V3), anchor = NA, processors = nc)

BrowseSeqs(lactob_type_16S_V1V3_aligned)

# V3-V4
lactob_type_16S_V3V4_aligned <-  DECIPHER::AlignSeqs(DNAStringSet(lb_type_16S_seqs_V3V4), anchor = NA, processors = nc)

BrowseSeqs(lactob_type_16S_V3V4_aligned)

# V4
lactob_type_16S_V4_aligned <-  DECIPHER::AlignSeqs(DNAStringSet(lb_type_16S_seqs_V4), anchor = NA, processors = nc)

BrowseSeqs(lactob_type_16S_V4_aligned)


# creating phylogenetic trees ---------------------------------------------

# WARNING CREATING PHYLOGENETIC TREES MAY TAKE A LOT OF TIME
# technically, I could have used a list with map instead of repeating the call
# to the function 4 times

# set this option to T if you plan to bootstrap the trees

boot_tree_option <- F

# phylogenetic tree, full sequences -----------------------------------------------
fitGTR_fullseq_lactobacillaceae <- cluster_seqs(alignment = lactob_type_16S_aligned, boot_option = boot_tree_option)

# phylogenetic tree, V1V3 ---------------------------------------------------------
fitGTR_V1V3_lactobacillaceae <- cluster_seqs(alignment = lactob_type_16S_V1V3_aligned, boot_option = boot_tree_option)

# phylogenetic tree, V3V4 ---------------------------------------------------------
fitGTR_V3V4_lactobacillaceae <- cluster_seqs(alignment = lactob_type_16S_V3V4_aligned, boot_option = boot_tree_option)

# phylogenetic tree, V4 ---------------------------------------------------------
fitGTR_V4_lactobacillaceae <- cluster_seqs(alignment = lactob_type_16S_V4_aligned, boot_option = boot_tree_option)

# bootstrapping the trees -------------------------------------------------

# this chunk optionally bootstraps the trees (***takes time***)
# depends on boot_tree_option, see above

# again, I am repeating the command 4 times, but the same can be done with map
# set mcore = F on Window


if(boot_tree_option){
  phylo_fullseq_lactobacillaceae <- boot_tree(fit = fitGTR_fullseq_lactobacillaceae)
  phylo_V1V3_lactobacillaceae <- boot_tree(fit = fitGTR_V1V3_lactobacillaceae)
  phylo_V3V4_lactobacillaceae <- boot_tree(fit = fitGTR_V3V4_lactobacillaceae)
  phylo_V4_lactobacillaceae <- boot_tree(fit = fitGTR_V4_lactobacillaceae)
}


# visualizing the 16S phylogenetic trees -----------------------------------------

# phylog. tree, full length 16S ---------------------------------------------
# optionally bootstrapping the tree
if(boot_tree_option){
  lbceae_ts_fulllength_tidytree <- as.treedata(phylo_fullseq_lactobacillaceae, boot = phylo_fullseq_lactobacillaceae$node.label)
} else {
  Lbceae_ts_full_length_16S_tree <- fitGTR_fullseq_lactobacillaceae$tree
  lbceae_ts_fulllength_tidytree <- as.treedata(Lbceae_ts_full_length_16S_tree, type = "ml")
 
}
# view the tibble
view(as_tibble(lbceae_ts_fulllength_tidytree))
# joining the metadata
class(seq_metadata)

lbceae_ts_fulllength_tidytree <- lbceae_ts_fulllength_tidytree %>% 
  left_join(seq_metadata, by = c("label" = "seq_abbr"))
view(as_tibble(lbceae_ts_fulllength_tidytree))

# the color scale
col_16S_group_df <- seq_metadata %>% select(Genus, group_color) %>% 
  arrange(Genus) %>% distinct()
# I am manually editing it
write_tsv(col_16S_group_df, "lactobacillaceae_colors.txt")
# modifiy and reopen if needed
# col_16S_group_df_ed <- read_tsv("lactobacillaceae_colors_edited.txt")
# otherwise
col_16S_group_df_ed <- col_16S_group_df
col_groups <-pull(col_16S_group_df_ed, group_color)
names(col_groups) <- pull(col_16S_group_df_ed, Genus)

# unrooted tree,  full length sequences ------------------------------------

if(boot_tree_option){
  lb16Sfl_treeplot_c <- ggtree(lbceae_ts_fulllength_tidytree, 
                                 aes(color = as.factor(Genus)),
                                 layout = "rectangular", 
                                 ladderize = T) + 
    geom_tiplab(aes(colour = Fermentation_type), size = 1) +
    geom_text(aes(label = bootstrap), hjust = 1, vjust = -0.4, size = 1, alpha = 0.5) +
    geom_treescale() +
    hexpand(0.2) +
    labs(title = "Lactobacillaceae, 16S RNA gene") +
    scale_colour_manual(values = col_groups) +
    guides(colour = "none") +
    theme_tree() +
    theme(plot.title = element_text(hjust = 0.5))
} else {
  lb16Sfl_treeplot_c <- ggtree(lbceae_ts_fulllength_tidytree, 
                               aes(color = as.factor(Genus)),
                               layout = "rectangular", 
                               ladderize = T) + 
    geom_tiplab(aes(colour = Fermentation_type), size = 1) +
    geom_treescale() +
    hexpand(0.2) +
    labs(title = "Lactobacillaceae, 16S RNA gene") +
    scale_colour_manual(values = col_groups) +
    guides(colour = "none") +
    theme_tree() +
    theme(plot.title = element_text(hjust = 0.5))
}

lb16Sfl_treeplot_c

# ggsave(lb16Sfl_treeplot_c, filename = str_c(taxo_group, "fl16S_c_unrooted.tiff", sep=""), dpi = 300)


# rooted tree, full length sequences --------------------------------------

# the root is Carnobacterium
outgroup_node <- lbceae_ts_fulllength_tidytree %>% as_tibble() %>%
  dplyr::filter(Genus == "Carnobacterium") %>% pull(node)
lbceae_ts_fulllength_tidytree_r <- treeio::root(phy = lbceae_ts_fulllength_tidytree, 
                                                outgroup = outgroup_node)

# b/w

lb16Sflr_treeplot <- ggtree(lbceae_ts_fulllength_tidytree_r, layout = "rectangular", ladderize = T) + 
  geom_tiplab(size = 1) +
  geom_treescale() +
  labs(title = "Lactobacillaceae, 16S RNA gene, rooted") +
  theme_tree() +
  theme(plot.title = element_text(hjust = 0.5))

# View(as.tibble(lbceae_ts_fulllength_tidytree_r))

if(boot_tree_option){
  lb16Sflr_treeplot_c <- ggtree(lbceae_ts_fulllength_tidytree_r, 
                                     aes(color = Genus), layout = "rectangular", 
                                     ladderize = T) + 
    geom_tiplab(size = 1) +
    geom_text(aes(label = bootstrap), hjust = 1, vjust = -0.4, size = 1, alpha = 0.5) +
    geom_treescale() +
    labs(title = "Lactobacillaceae, 16S RNA gene, rooted") +
    scale_color_manual(values = col_groups) +
    guides(color = "none") +
    theme_tree() +
    theme(plot.title = element_text(hjust = 0.5))


} else {
  
  lb16Sflr_treeplot_c <- ggtree(lbceae_ts_fulllength_tidytree_r, 
                                aes(color = Genus), layout = "rectangular", 
                                ladderize = T) + 
    geom_tiplab(size = 1) +
    geom_treescale() +
    labs(title = "Lactobacillaceae, 16S RNA gene, rooted") +
    scale_color_manual(values = col_groups) +
    guides(color = "none") +
    theme_tree() +
    theme(plot.title = element_text(hjust = 0.5))

}

lb16Sflr_treeplot_c
# the following is supplementary figure 1 in the manuscript
# ggsave(lb16Sflr_treeplot, filename = str_c(taxo_group, "fl16S_r.tiff", sep=""), dpi = 300)

# ggsave(lb16Sflr_treeplot_c, filename = str_c(taxo_group, "fl16S_r_color.pdf", sep=""),
#       height = 9, width = 8)


# plotting tree V1-V3 --------------------------------------------------------
# if you want to add bootstrap values copy the code for full sequences
if(boot_tree_option){
  lbceae_ts_V1V3_tidytree <- as.treedata(phylo_V1V3_lactobacillaceae, boot = phylo_V1V3_lactobacillaceae$node.label)
} else {
  lbceae_ts_V1V3_16S_tree <- fitGTR_V1V3_lactobacillaceae$tree
  lbceae_ts_V1V3_tidytree <- as.treedata(lbceae_ts_V1V3_16S_tree, type = "ml")
}

# joining the metadata
lbceae_ts_V1V3_tidytree <- lbceae_ts_V1V3_tidytree %>% left_join(seq_metadata, by = c("label" = "seq_abbr"))
# view(as_tibble(lbceae_ts_V1V3_tidytree))
# reroot and plot
outgroup_node_V1V3 <- lbceae_ts_V1V3_tidytree %>% as_tibble() %>%
  dplyr::filter(Genus == "Carnobacterium") %>% pull(node)
lbceae_ts_V1V3_tidytree_r <- treeio::root(phy = lbceae_ts_V1V3_tidytree, 
                                          outgroup = outgroup_node_V1V3)
# view(as_tibble(lbceae_ts_V1V3_tidytree_r))

lb16SV1V3r_treeplot_c <- ggtree(lbceae_ts_V1V3_tidytree_r, 
                                aes(color = Genus),
                                layout = "rectangular", ladderize = T) + 
  geom_tiplab(mapping = aes(colour = Fermentation_type), size = 1) +
  geom_treescale() +
  labs(title = "Lactobacillaceae, 16S RNA gene, V1V3, rooted") +
  scale_color_manual(values = col_groups) +
  guides(color = "none") +
  theme_tree() +
  theme(plot.title = element_text(hjust = 0.5))
lb16SV1V3r_treeplot_c

# ggsave(lb16SV1V3r_treeplot_c, filename = str_c(taxo_group, "V1V316S_r_c.tiff", sep=""), dpi = 300,
#       height = 9, width = 8)


# plotting tree V3-V4 --------------------------------------------------------
# if you want to add bootstrap values copy the code for full sequences
if(boot_tree_option){
  lbceae_ts_V3V4_tidytree <- as.treedata(phylo_V3V4_lactobacillaceae, boot = phylo_V3V4_lactobacillaceae$node.label)
} else {
  lbceae_ts_V3V4_16S_tree <- fitGTR_V3V4_lactobacillaceae$tree
  lbceae_ts_V3V4_tidytree <- as.treedata(lbceae_ts_V3V4_16S_tree, type = "ml")
}

# view the tibble
# view(as_tibble(lbceae_ts_V3V4_tidytree))
# joining the metadata
lbceae_ts_V3V4_tidytree <- left_join(lbceae_ts_V3V4_tidytree, 
                                     seq_metadata, by = c("label" = "seq_abbr"))
# view(as_tibble(lbceae_ts_V3V4_tidytree))
# reroot and plot
outgroup_node_V3V4 <- lbceae_ts_V3V4_tidytree %>% as_tibble() %>%
  dplyr::filter(Genus == "Carnobacterium") %>% pull(node)
lbceae_ts_V3V4_tidytree_r <- treeio::root(phy = lbceae_ts_V3V4_tidytree, 
                                          outgroup = outgroup_node_V3V4)
# view(as_tibble(lbceae_ts_V3V4_tidytree_r))

lb16SV3V4r_treeplot_c <- ggtree(lbceae_ts_V3V4_tidytree_r, 
                                mapping = aes(colour = Genus),
                                layout = "rectangular", ladderize = T) + 
  geom_tiplab(mapping = aes(colour = Fermentation_type), size = 1) +
  geom_treescale() +
  labs(title = "Lactobacillaceae, 16S RNA gene, V3V4, rooted") +
  scale_color_manual(values = col_groups) +
  guides(color = "none") +
  theme_tree() +
  theme(plot.title = element_text(hjust = 0.5))
lb16SV3V4r_treeplot_c

# ggsave(lb16SV3V4r_treeplot_c, filename = str_c(taxo_group, "V3V416S_r_c.tiff", sep=""), 
#        dpi = 300, height = 9, width = 8)


# plotting tree  V4 --------------------------------------------------------
# if you want to add bootstrap values copy the code for full sequences

if(boot_tree_option){
  lbceae_ts_V4_tidytree <- as.treedata(phylo_V4_lactobacillaceae, boot = phylo_V4_lactobacillaceae$node.label)
} else {
  lbceae_ts_V4_16S_tree <- fitGTR_V4_lactobacillaceae$tree
  lbceae_ts_V4_tidytree <- as.treedata(lbceae_ts_V4_16S_tree, type = "ml")
}

# view the tibble
# view(as_tibble(lbceae_ts_V4_tidytree))
# joining the metadata
lbceae_ts_V4_tidytree <- left_join(lbceae_ts_V4_tidytree, 
                                   seq_metadata, by = c("label"="seq_abbr"))
# view(as_tibble(lbceae_ts_V4_tidytree))
# reroot and plot
outgroup_node_V4 <- lbceae_ts_V4_tidytree %>% as_tibble() %>%
  dplyr::filter(Genus == "Carnobacterium") %>% pull(node)
lbceae_ts_V4_tidytree_r <- treeio::root(phy = lbceae_ts_V4_tidytree, 
                                        outgroup = outgroup_node_V4)
# view(as_tibble(lbceae_ts_V4_tidytree_r))

lb16SV4r_treeplot_c <- ggtree(lbceae_ts_V4_tidytree_r, 
                              mapping = aes(colour = Genus),
                              layout = "rectangular", ladderize = T) + 
  geom_tiplab(aes(colour = Fermentation_type), size = 1) +
  geom_treescale() +
  labs(title = "Lactobacillaceae, 16S RNA gene, V4, rooted") +
  scale_color_manual(values = col_groups) +
  guides(colour = "none") +
  theme_tree() +
  theme(plot.title = element_text(hjust = 0.5))
lb16SV4r_treeplot_c

# ggsave(lb16SV4r_treeplot_c, filename = str_c(taxo_group, "V416S_r_c.tiff", sep=""), dpi = 300,
#       height = 9, width = 8)


# simplifying the tree ----------------------------------------------------

# collapse clades

# get the parent node for genus Lactobacillus
Lactobacillus_tips <- lbceae_ts_fulllength_tidytree_r %>% as_tibble() %>%
  dplyr::filter(Genus == "Lactobacillus") %>% pull(node)
MRCA_lactobacillus <- MRCA(lbceae_ts_fulllength_tidytree_r, Lactobacillus_tips)
# this probably needs to be done with a function using map or a loop
# once that is done I can annotate the corresponding node with the appropriate label
# may be using geom_cladelab
# it is interesting to notices that geom_cladelab takes a secondary data frame

# can also be used to create a tree for a clade
viewClade(lb16Sflr_treeplot_c, MRCA_lactobacillus)

# let's see if I can get a df with MRCAs using a loop
# let's exploit the species per genus df, which already contains colors
# I need to exclude the outgroup
species_per_genus <- seq_metadata %>% group_by(Genus) %>% summarize(n_species = n())
MRCA_metadata <- seq_metadata %>% 
  dplyr::select(Genus, Fermentation_type, group_16S, group_color) %>%
  distinct() %>% 
  left_join(select(species_per_genus, Genus, n_species))

# need a workaround because of Leuc. fallax, and Weissella muntiaci 
MRCAs_16S_fl <- vector(mode = "list", length = nrow(MRCA_metadata))
for (i in seq_along(MRCA_metadata$Genus)){
  tips_temp <- lbceae_ts_fulllength_tidytree_r %>% as_tibble() %>%
    dplyr::filter(Genus == MRCA_metadata$Genus[[i]]) 
  if(!any(tips_temp$Genus %in% c("Leuconostoc", "Lentilactobacillus", "Weissella")) ) {
    # this will handle correctly all genera except the problematic ones
    tips_temp_v <-  pull(tips_temp, node) 
  } else {
    tips_temp_v <-  tips_temp %>% 
      dplyr::filter(!(Species %in% c("fallax","muntiaci","oryzae","kribbianus", "kosonis", "senioris", "curieae"))) %>% 
      pull(node) 
    # this handle special cases 
  }
  MRCAs_16S_fl[[i]] <- MRCA(lbceae_ts_fulllength_tidytree_r, tips_temp_v)
  names(MRCAs_16S_fl)[i] <- MRCA_metadata$Genus[[i]]
}
rm(tips_temp, tips_temp_v)
MRCAs_16S_fl_df <- as.data.frame(do.call(rbind, MRCAs_16S_fl)) %>% rownames_to_column(var ="Genus")
MRCA_metadata <- left_join(MRCA_metadata, MRCAs_16S_fl_df) %>%
  dplyr::rename(MRCAnc = V1) %>%
  mutate(is_clade = (n_species>1))


# let's try geom_cladelab
df <- dplyr::rename(MRCA_metadata, genus = Genus, Ftype = Fermentation_type)

lb16Sflr_treeplot_2 <- ggtree(lbceae_ts_fulllength_tidytree_r, layout = "rectangular", ladderize = T) + 
  geom_treescale() +
  labs(title = "Lactobacillaceae, 16S RNA gene, rooted") +
  theme_tree() +
  theme(plot.title = element_text(hjust = 0.5))


cladetree_16Sfl_r <- lb16Sflr_treeplot_2 + 
  geom_cladelab(data = df,
                mapping = aes(node = MRCAnc,
                              label = genus,
                              color = Ftype), 
                size = 0.25, fontface = 4) +
  guides(color = "none") +
  hexpand(.1) +
  vexpand(.05, direction = -1)
cladetree_16Sfl_r  

# this is Figure 1 in the manuscript
ggplot2::ggsave(cladetree_16Sfl_r, filename = "cladelabeled_fl16S.pdf",
                 width = 12, height = 18)
ggplot2::ggsave(cladetree_16Sfl_r, filename = "cladelabeled_fl16S.tiff",
                 width = 12, height = 18, dpi = 300)

# the tree and the alignment
msaplot(cladetree_16Sfl_r, "lactob_type_16S_aligned.fasta")


# kmer analysis ------------------------------------------------------------

# using as.DNAbin.character()
# let's start by creating a DNAbin object from the unaligned matrix
lb_type_16S_type_fl <- ape::read.FASTA("lactob_type_16S.fasta")


# lb_type_16S_type_fl

# creating the distance matrix using pentamers (can easily increase, it is very fast)

lb_type_16S_type_fl.kdist <- kdistance(lb_type_16S_type_fl, k = 5) # very quick
print(as.matrix(lb_type_16S_type_fl.kdist)[1:10, 1:10], digits = 3)

dim(as.matrix(lb_type_16S_type_fl.kdist))

otus_lb_16S_fl_097_km5 <- otu(lb_type_16S_type_fl, k = 5, threshold = 0.97, method = "farthest", nstart = 20)
# this takes a few seconds
length(unique(otus_lb_16S_fl_097_km5))

# only 184, probably because closely related species are merged in a single otu
# I am using the 0.9865 threshold for species suggested by Qiao, N., Wittouck, S., 
# Mattarelli, P., Zheng, J., Lebeer, S., Felis, G.E., Gänzle, M.G., 2022. 
# After the storm—Perspectives on the taxonomy of Lactobacillaceae. Jds Commun 3, 
# 222–227. https://doi.org/10.3168/jdsc.2021-0183
otus_lb_16S_fl_099_km5 <- otu(lb_type_16S_type_fl, k = 5, threshold = 0.9865, method = "farthest", nstart = 20)
# this takes a few seconds
length(unique(otus_lb_16S_fl_099_km5))
# 285
otus_lb_16S_fl_099_km5[1:10]
# convert to data frames, and removing a final asterisk (technically can be merged with sequence metadata)
otus_lb_16S_fl_099_km5 <- tibble::enframe(otus_lb_16S_fl_099_km5) %>%
  rename(seq_abbr = name, OTU_99_km5 = value) %>%
  mutate(seq_abbr = str_remove(seq_abbr, "\\*"))

# genera with a larger number of species gathered in a single OTUs
genera_w_cl_rel_species <- left_join(otus_lb_16S_fl_099_km5, seq_metadata) %>%
  group_by(Genus, OTU_99_km5) %>%
  summarize (n = n()) %>%
  dplyr::arrange(Genus, desc(n))
write_tsv(genera_w_cl_rel_species, file = "genera_w_cl_rel_species_fl.txt")

# now, this is probably going to be useful when working with reference sequences from SILVA


# MDS with vegan using kmer distance matrix -------------------------------

# this and metaMDS fail MDS_lb_16S_fl_km5 <- MASS::isoMDS(lb_type_16S_type_fl.kdist, maxit = 50)
# this is a PCOA
MDS_lb_16S_fl_km5 <- cmdscale(lb_type_16S_type_fl.kdist, list. =T)
MDS_lb_16S_fl_km5_coord <- as.data.frame(MDS_lb_16S_fl_km5$points) %>%
  rownames_to_column(var = "seq_abbr") %>%
  rename(dim1 = V1, dim2 = V2)
# join info from seq_metadata
MDS_lb_16S_fl_km5_coord <- left_join(MDS_lb_16S_fl_km5_coord,
                                     dplyr::select(seq_metadata, seq_abbr,
                                                   s_label:Accn_n, 
                                                   ferm_type = Fermentation_type,
                                                   group_16S, group_color))
# creating a df with centroids for genera
MDS_lb_16S_fl_km5_coord_centr <- MDS_lb_16S_fl_km5_coord %>%
  group_by(Genus, group_color, ferm_type) %>%
  summarize(mdim1 = mean(dim1),
            mdim2 =  mean(dim2)) 

PCoA_plot_16S_fl_km5 <- MDS_lb_16S_fl_km5_coord %>%
  ggplot(mapping = aes(x = dim1, y = dim2, color = Genus)) +
  geom_point() + 
  scale_color_manual(values = col_groups) +
  labs(x = "dim(1)",
       y = "dim(2)") +
  theme_bw()

PCoA_plot_16S_fl_km5
# the display is not very useful except in showing there is a lot of overlap

# kmer V1-V3
# let's start by creating a DNAbin object from the unaligned matrix

lb_type_16S_type_V1V3 <- insect::char2dna(lb_type_16S_seqs_V1V3)

lb_type_16S_type_V1V3
# creating the distance matrix using pentamers (can easily increase, it is very fast)

lb_type_16S_type_V1V3.kdist <- kdistance(lb_type_16S_type_V1V3, k = 5) # very quick
# print(as.matrix(lb_type_16S_type_V1V3.kdist)[1:10, 1:10], digits = 3)

# dim(as.matrix(lb_type_16S_type_V1V3.kdist))

otus_lb_16S_V1V3_097_km5 <- otu(lb_type_16S_type_V1V3, k = 5, threshold = 0.97, method = "farthest", nstart = 20)
# this takes a few seconds
length(unique(otus_lb_16S_V1V3_097_km5))
# only 254, probably because several closely related species are merged in a single otu
otus_lb_16S_V1V3_099_km5 <- otu(lb_type_16S_type_V1V3, k = 5, threshold = 0.9865, method = "farthest", nstart = 20)
# this takes a few seconds
length(unique(otus_lb_16S_V1V3_099_km5))
# 320
otus_lb_16S_V1V3_099_km5[1:10]
# convert to data frames, and removing a final asterics (technically can be merged with sequence metadata)
otus_lb_16S_V1V3_099_km5 <- tibble::enframe(otus_lb_16S_V1V3_099_km5) %>%
  rename(seq_abbr = name, OTU_99_km5 = value) %>%
  mutate(seq_abbr = str_remove(seq_abbr, "\\*"))

otus_lb_16S_V1V3_099_km5 <- otus_lb_16S_V1V3_099_km5 %>%
  left_join(dplyr::select(seq_metadata,
                          seq_abbr,
                          s_label:Accn_n, 
                          ferm_type = Fermentation_type,
                          group_16S, group_color
  ))

# write_tsv(otus_lb_16S_V1V3_099_km5, file = "otus_lb_16S_V1V3_099_km5.txt")
# genera with a larger number of species gathered in a single OTUs
genera_w_cl_rel_species_V1V3 <- left_join(otus_lb_16S_V1V3_099_km5, seq_metadata) %>%
  group_by(Genus, OTU_99_km5) %>%
  summarize (n = n()) %>%
  dplyr::arrange(Genus, desc(n))
write_tsv(genera_w_cl_rel_species_V1V3, file = "genera_w_cl_rel_species_V1V3.txt")

# kmer V3-V4
# let's start by creating a DNAbin object from the unaligned matrix
lb_type_16S_type_V3V4 <- insect::char2dna(lb_type_16S_seqs_V3V4)
lb_type_16S_type_V3V4
# creating the distance matrix using pentamers (can easily increase, it is very fast)

lb_type_16S_type_V3V4.kdist <- kdistance(lb_type_16S_type_V3V4, k = 5) # very quick
# print(as.matrix(lb_type_16S_type_V3V4.kdist)[1:10, 1:10], digits = 3)

# dim(as.matrix(lb_type_16S_type_V3V4.kdist))

otus_lb_16S_V3V4_097_km5 <- otu(lb_type_16S_type_V3V4, k = 5, threshold = 0.97, method = "farthest", nstart = 20)
# this takes a few seconds
length(unique(otus_lb_16S_V3V4_097_km5))
# only 149, probably because several closely related species are merged in a single otu
otus_lb_16S_V3V4_099_km5 <- otu(lb_type_16S_type_V3V4, k = 5, threshold = 0.9865, method = "farthest", nstart = 20)
# this takes a few seconds
length(unique(otus_lb_16S_V3V4_099_km5))
# 211
otus_lb_16S_V3V4_099_km5[1:10]
# convert to data frames, and removing a final asterisc (technically can be merged with sequence metadata)
otus_lb_16S_V3V4_099_km5 <- tibble::enframe(otus_lb_16S_V3V4_099_km5) %>%
  rename(seq_abbr = name, OTU_99_km5 = value) %>%
  mutate(seq_abbr = str_remove(seq_abbr, "\\*"))

otus_lb_16S_V3V4_099_km5 <- otus_lb_16S_V3V4_099_km5 %>%
  left_join(dplyr::select(seq_metadata,
                          seq_abbr,
                          s_label:Accn_n, 
                          ferm_type = Fermentation_type,
                          group_16S, group_color
  ))

# write_tsv(otus_lb_16S_V3V4_099_km5, file = "otus_lb_16S_V3V4_099_km5.txt")
genera_w_cl_rel_species_V3V4 <- left_join(otus_lb_16S_V3V4_099_km5, seq_metadata) %>%
  group_by(Genus, OTU_99_km5) %>%
  summarize (n = n()) %>%
  dplyr::arrange(Genus, desc(n))
write_tsv(genera_w_cl_rel_species_V3V4, file = "genera_w_cl_rel_species_V3V4.txt")

# kmer V4
# let's start by creating a DNAbin object from the unaligned matrix
lb_type_16S_type_V4 <- insect::char2dna(lb_type_16S_seqs_V4)
lb_type_16S_type_V4
# creating the distance matrix using pentamers (can easily increase, it is very fast)

lb_type_16S_type_V4.kdist <- kdistance(lb_type_16S_type_V4, k = 5) # very quick
# print(as.matrix(lb_type_16S_type_V4.kdist)[1:10, 1:10], digits = 3)

# dim(as.matrix(lb_type_16S_type_V4.kdist))

otus_lb_16S_V4_097_km5 <- otu(lb_type_16S_type_V4, k = 5, threshold = 0.97, method = "farthest", nstart = 20)
# this takes a few seconds
length(unique(otus_lb_16S_V4_097_km5)) 
# 114
# only 114, probably because several closely related species are merged in a single otu
otus_lb_16S_V4_099_km5 <- otu(lb_type_16S_type_V4, k = 5, threshold = 0.9865, method = "farthest", nstart = 20)
# this takes a few seconds
length(unique(otus_lb_16S_V4_099_km5))
# 174
otus_lb_16S_V4_099_km5[1:10]
# convert to data frames, and removing a final asterik (technically can be merged with sequence metadata)
otus_lb_16S_V4_099_km5 <- tibble::enframe(otus_lb_16S_V4_099_km5) %>%
  rename(seq_abbr = name, OTU_99_km5 = value) %>%
  mutate(seq_abbr = str_remove(seq_abbr, "\\*"))

otus_lb_16S_V4_099_km5 <- otus_lb_16S_V4_099_km5 %>%
  left_join(dplyr::select(seq_metadata,
                          seq_abbr,
                          s_label:Accn_n, 
                          ferm_type = Fermentation_type,
                          group_16S, group_color
  ))
# write_tsv(otus_lb_16S_V4_099_km5, file = "otus_lb_16S_V4_099_km5.txt")
# genera with a larger number of species gathered in a single OTUs
genera_w_cl_rel_species_V4 <- left_join(otus_lb_16S_V4_099_km5, seq_metadata) %>%
  group_by(Genus, OTU_99_km5) %>%
  summarize (n = n()) %>%
  dplyr::arrange(Genus, desc(n))
write_tsv(genera_w_cl_rel_species_V4, file = "genera_w_cl_rel_species_V4.txt")





# I am creating a function which:
# a. takes as an input a named character vector with sequences (argument seq_char)
# b. converts it to DNA.bin
# c. uses functions of package kmer with kmer_n length (default = 5) to
#    1. calculate kmer distance matrix
#    2. return OTU membership vector at two similarity levels (sim_lvl_1 = 0.97, sim_lvl_2 = 0.99)
#    3. optionally perform MDS and return results and coordinates (which otherwise are null)
# results ar returned as a list

# I am not doing any effort for at error trapping (not even checking if required packages
# have been loaded: anyway kmer is needed and also ape, as a dependency, but it is loaded
# by phylotools)

kmer_analysis <- function(seq_char, kmer_n = 5, 
                          sim_lvl_1 = 0.97, sim_lvl_2 = 0.9865,
                          perform_MDS = T){
  # create DNA_bin object
  seq_DNAbin <- insect::char2dna(seq_char)
  # creating the distance matrix
  kdist_mat <- kdistance(seq_DNAbin, k = kmer_n) 
  
  otus_lvl1 <- otu(seq_DNAbin, k = 5, threshold = sim_lvl_1, method = "farthest", nstart = 20)
  
  otus_lvl1_df <- tibble::enframe(otus_lvl1) %>%
    rename(seq_abbr = name, OTU = value) %>%
    mutate(seq_abbr = str_remove(seq_abbr, "\\*"))
  
  otus_lvl2 <- otu(seq_DNAbin, k = 5, threshold = sim_lvl_2, method = "farthest", nstart = 20)
  
  otus_lvl2_df <- tibble::enframe(otus_lvl2) %>%
    rename(seq_abbr = name, OTU = value) %>%
    mutate(seq_abbr = str_remove(seq_abbr, "\\*"))
  
  if(perform_MDS) {
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


# Citations ---------------------------------------------------------------

map(c(.cran_packages, .bioc_packages), citation)

# Credits and copyright ---------------------------------------------------


# Assume that this is overall under MIT licence

# Copyright 2022 Eugenio Parente
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
