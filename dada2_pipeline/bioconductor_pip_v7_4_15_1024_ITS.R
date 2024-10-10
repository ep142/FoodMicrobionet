# DADA2/Bioconductor pipeline for ITS, modified, v7.4.15, 4/10/24

#  Description & instructions ---------------------------------------------

# This script is designed to process sequences for small studies (<60-80 
# samples, but even 100  can be processed depending on RAM) of amplicon targeted 
# metagenomic data for fungi using the DADA2 pipeline 
# https://benjjneb.github.io/dada2/tutorial.html
# with modifications for ITS analysis 
# https://benjjneb.github.io/dada2/ITS_workflow.html
# and to carry out further processing needed to ready the output for import
# into FoodMicrobionet (https://github.com/ep142/FoodMicrobionet)

# This version of the script includes options for both data downloaded from NCBI
# SRA and data obtained from Novogene UK Ltd. There are also options for
# single end/paired end data sets obtained with Illumina or 454 or Ion Torrent
# In addition it includes a workaround for binned quality fastq described here:
# https://github.com/benjjneb/dada2/issues/791 and a workaround for  handling
# ITS sequences which do not merge properly due to excess length described here:
# https://github.com/benjjneb/dada2/issues/537
# in this version UNITE general release for fungi with singletons is used as a taxonomic reference
# Finally, the script assembles an object which can be used in the future for
# redoing the taxonomic assignment in an automated way

# to use this script

# 1. copy the script to a new folder and create a RStudio project; the directory 
#    containing the project folder must also contain the bioconductor_pip_ITS_functions.R
#    file and the primer_pairs_fungi.txt file (see below)
# 2. create a new folder called data
# 3. if you are using data downloaded from NCBI SRA
#    3.1 inside the folder data create three folders: fastq, filtered, metadata
#    3.2 download from SRA the accession list and the study metadata for the
#        any study you want to process (must be ITS) 
#    3.3 using the sratoolkit download the fastq files from SRA and put them in the fastq folder
#    3.4 put the metadata in the metadata folder
# 4. if you are using data delivered by Novogene
#    4.1 put the folder you received from Novogene in the project folder
#    4.2 create the filtered and metadata folders within the data folder
#    4.3 put the metadata in the metadata folder; the first column must match 
#        sample names (as received from Novogene, these are also the names of 
#        the folders containing the sequences)
# 5. download the most recent taxonomy reference files from 
#    https://benjjneb.github.io/dada2/training.html
#    and put them in a folder called taxdb at the same level of the directory containing your project
# 6. an optional table with commonly used primer pairs (primer_pairs_fungi.txt) 
#    can be downloaded from our repository and is used by a function to 
#    check primers; it should be placed in the folder containing the folder of 
#.   your project and the taxonomic reference folder.
# 7. make sure you have all the information you need (primers, platform, region, etc.)

# If you want to use cutadapt to remove primers install it using miniconda
# A good tutorial is here https://astrobiomike.github.io/unix/conda-intro
# You can easily adapt the script to your own data
# If you want to give it a try, download from SRA the study used in section study metadata
# Run the script one section at a time and check the output and warnings, if any

# Install/load packages ---------------------------------------------------

.cran_packages <- c("tidyverse", "parallel", "beepr", "tictoc", "reticulate", 
                    "logr", "R.utils", "crayon")
.bioc_packages <- c("BiocManager","dada2", "phyloseq", "DECIPHER", "phangorn", 
                    "BiocStyle", "ShortRead")

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

# the following command detects the number of cores on UNIX/MacOS
nc <- parallel::detectCores(logical = F) # to detect physical cores in MacOS
set.seed(100)
# mainly for reproducibility reasons (this will be saved with the workspace)
r_version <- R.Version() # in the future may be check that the version running is compatible
# with the script (it only matters if somebody other than myself is using the script)
session_info <- sessionInfo()

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

# standard options for handling primers: you have to change this manually 
# if you want to do otherwise
check_primers <- "auto" # "visually" if you want to do visually, "auto" otherwise
use_cutadapt <- T # use cutadapt to remove primers, takes time and space, use wisely
use_primer_table <- T # uses a primer table to check if the primer sequences are OK
# put the location of your primer table inside this if
# hint: if the file is two levels up just duplicate ".." below
if(use_primer_table){
  primer_table_path <- file.path("..", "primer_pairs_fungi.txt")
}

# functions ---------------------------------------------------------------

# this will load three functions which are going to be used in the pipeline
# this file should be in the folder containing your project directory
# otherwise you need to change the path
# hint: if the source files is two levels up just duplicate ".." below
source(file.path("..","bioconductor_pip_ITS_functions.R"))

# study metadata ----------------------------------------------------------

# change this as appropriate
data_type <- "sra" 
# alternatives are "sra", for data downloaded from sra
# "novogene_raw", for data obtained from novogene, with primer not removed
# "novogene_rawa", for data obtained from novogene, with primer  and adapter not removed
# "novogene_clean", for data obtained from novogene, data with primer removed, 
#  after merging and chimera removal

# creating information for the study and sample data frames
# when using data other than those downloaded from SRA replace
# the accession number with whichever identifier you want to use

Study <- "SRP108314" # PRJNA388367
target <- "ITS region and 16S RNA gene" # or ITS region and 16S RNA gene (or 16S RNA)
region <- "ITS2 and V3-V4" # or ITSx and Vx
seq_accn <- Study
DOI <- "10.1016/j.micres.2017.09.004"

# information on the platform and arrangement

platform <- "Illumina_miseq" # (or set to "Illumina_miseq", "Illumina_novaseq", "Illumina_HiSeq", "Illumina_iSeq" or "Ion_Torrent" or "F454")
# 
paired_end <- T # set to true for paired end, false for single end or in the case of clean, merged seqs
if (!paired_end) overlapping <- T # needed to run species assignment for SILVA

# Primers --------------------------------------------
# all in standard 5'-3' orientation
# ITS1:
# ITS1-F_KYO2 (18S SSU 1733–1753) TAGAGGAAGTAAAAGTCGTAA 21 bp and ITS2_KYO2 (5.8 2046–2029) CTHGGTCATTTAGAGGAASTAA 22 bp
# (Toju et al., 2012)
# ITS1FI2 GAACCWGCGGARGGATCA (18 bp) and 5.8S CGCTGCGTTCTTCATCG (17 bp)
# BITS ACCTGCGGARGGATCA (18 bp) and B58S3 GAGATCCRTTGYTRAAAGTT (20 bp)  (Bokulich and Mills, 2013)
# ITS1 TCCGTAGGTGAACCTGCGG (19 bp) and ITS4 TCCGTAGGTGAACCTGCGG (19 bp)
# ITS5F GGAAGTAAAAGTCGTAACAAGG (22 bp)	ITS1R	GCTGCGTTCTTCATCGATGC (20 bp)
# ITS1F TTGGTCATTTAGAGGAAGTAA (21 bp) ITS2 GCTGCGTTCTTCATCGATGC (20 bp)
# ITS1Fv2 CTTGGTCATTTAGAGGAAGTAA (22 bp) ITS1R GCTGCGTTCTTCATCGATGC (20 bp) (White et al., 1990)

# ITS2: 
# ITS3F GCATCGATGAAGAACGCAGC ITS4R TCCTCCGCTTATTGATATGC 
# (White TJ, Bruns TD, Lee SB, Taylor JW (1990) Amplification and direct sequencing of fungal ribosomal RNA genes for phylogenetics. In: Innis MA,Gelfand DH, Sninsky JJ, White TJ, editors. PCR protocols: a guide to methodsand applications. United States: Academic Press. pp. 315–322)
# ITS86F GTGAATCATCGAATCTTTGAA (21 bp) and ITS4 TCCTCCGCTTATTGATATGC (20 bp)
# ITS3f GCATCGATGAAGAACGCAGC (20 bp) and ITS4-KYO1 TCCTCCGCTTWTTGWTWTGC (20 bp) Toju et al., 2012
# F2045 GCATCGATGAAGAACGCAGC (20 bp) and R2390 TCCTCCGCTTATTGATATGC (20 bp)
# ITS3-KYO2 GATGAAGAACGYAGYRAA (18 bp)) and ITS4 TCCTCCGCTTATTGATATGC (20 bp)

# ITS1 5.8S and ITS2:
# ITS1 TCCGTAGGTGAACCTGCGG (19 bp) and ITS4 TCCTCCGCTTATTGATATGC (20 bp)

# expected amplicon length is variable
primer_f <- "ITS3-KYO2" # 18 bp
primer_r <- "ITS4"  # 20 bp

# NOTE: be extra careful in indicating primer sequences because this will affect
# primer detection and primer removal by cutadapt

FWD <- "GATGAAGAACGYAGYRAA"  ## CHANGE THIS to your forward primer sequence
REV <- "TCCTCCGCTTATTGATATGC" ## CHANGE THIS to your reverse primer sequence

target1 <- "ITS_DNA"
target2 <- region

# if you have set correctly your options this should not return any warning
# but will return details on the primer table in the console
if(use_primer_table) double_check_primers(
  primerf_name = primer_f, primerr_name = primer_r,
  primerf_seq = FWD, primerr_seq = REV
)

# ad hoc: in this case the region covered by the primers is quite long
# be careful with merging and if it fails set overlapping <- F

(FWD.orients <- allOrients(FWD))
(REV.orients <- allOrients(REV))

# create sub-directory ----------------------------------------------------------
# sub-directory data must already be in the wd
if(data_type == "sra"){
  fastq_path <- file.path("data", "fastq")
  file_list <- list.files(fastq_path)
  fns <- sort(list.files(fastq_path, full.names = TRUE))
} else{
  # this assumes that the only alternative is data from novogene
  top_level_dir <- list.dirs(recursive = F)[str_detect(list.dirs(recursive = F), "result")]
  data_dirs <- list.dirs(file.path(top_level_dir), recursive = F)
  if(str_detect(data_type, "clean")){
    paired_end<-F
    fastq_dir <- data_dirs[str_detect(data_dirs,"Clean")]
  } else {
    paired_end<-T
    fastq_dir <- data_dirs[str_detect(data_dirs,"Raw")]
  }
  fastq_dirs <- list.dirs(fastq_dir, recursive = F)
  file_list <- vector("list",length = length(fastq_dirs))
  for(i in seq_along(fastq_dirs)){
    file_list[[i]]<-list.files(fastq_dirs[i])
    if(str_detect(data_type, "clean")){
      file_index <- which(str_detect(file_list[[i]], "effective.fastq.gz"))
    } else {
      if(str_detect(data_type, "rawa")){
        file_index <- which(str_detect(file_list[[i]], "raw"))
      } else {
        file_index <- which(!str_detect(file_list[[i]], "raw") & !str_detect(file_list[[i]], "Frags"))
      }
    }
    file_list[[i]]<-file.path(
      fastq_dirs[i],
      file_list[[i]][file_index]
    )
  }
  fns <- unlist(file_list)
  file_list <- basename(unlist(file_list))
}

filt_path <- file.path("data", "filtered")

# if(!file_test("-d", fastq_path)) {
#  dir.create(fastq_path)
# only needed if you want to create the directory and move the data from somewhere else

if(paired_end){
  fnFs <- fns[grepl("_1", fns)]
  fnRs <- fns[grepl("_2", fns)]
} else {
  fnFs <- fns
}

# separates the name of the samples from the names of the files; must be adapted
# for different datasets
# the separator for names is "." for 454 (this only removes extension) and generally "_" for Illumina
sample.names <- dplyr::case_when(
  all(str_detect(basename(fnFs),"_")) ~ sapply(strsplit(basename(fnFs), "_", fixed = T), `[`, 1),
  all(str_detect(basename(fnFs),".")) ~ sapply(strsplit(basename(fnFs), ".", fixed = T), `[`, 1),
)
# ad hoc
if(data_type == "novogene_rawa"){
  sample.names <- str_remove(sample.names, "\\.raw")
}

# early check for sequences and metadata --------------------------------

# path to the metadata file (needs to be adapted)
metadata_path <- file.path("data", "metadata", "SraRunTable.txt") 
# alternatives when using data downloads from NCBI SRA are:
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
    dplyr::slice(1:18) |>
    dplyr::select(Run)
  write_tsv(acc_list_bacteria, file = "acc_list_bacteria.txt")
  acc_list_fungi <- samdf |>
    arrange(Run) |>
    dplyr::slice(19:35) |>
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
  if(all(sample.names %in% samdf$Run)){
    cat("\nsamples in fastq files match samples in metadata\n")
  } else {
    cat("\nWARNING samples in fastq files DO NOT match samples in metadata\n")
  }
} else{
  # ad hoc for your data
  if(all(sample.names %in% samdf$Library_Name)){
    cat("\nsamples in fastq files match samples in metadata\n")
  } else {
    cat("\nWARNING samples in fastq files DO NOT match samples in metadata\n")
  }
}


# handling primers --------------------------------------------------------
# change this if you do not want to use standard options
# check for the occurrence of primers visually or automatically
# use visually if you want to check visually  or auto as an alternative, it is much more efficient
# check_primers <- "auto" 

# sequences will sampled and shown here as a double check

# check carefully if primer occur in the right position: 
# the forward primer  should appear at the beginning of forward sequences 
# the reverse primer should occur at the beginning of reverse sequences (if any)
# extra nt may occasionally be present
# if you have chosen to check primer visually take note now of their occurrence 
# and the length of the forward and reverse sequence to trim at the 5' end

if(keep_time) tic("\nReading sequences")
# check for occurrence of primers and adapters on a sample of forward and
# reverse sequences
# get a sample of 6 sequences
sampleFs <- if(length(fnFs)>=6) sample(fnFs,6) else fnFs[1:length(fnFs)]
# works with .fastq and fastq.gz
myFwsample <- ShortRead::readFastq(sampleFs)
print(head(sread(myFwsample),10))
print(tail(sread(myFwsample),10))
ave_seq_length_f <- round(mean(sread(myFwsample)@ranges@width))
ave_seq_length <- ave_seq_length_f

cat("\naverage sequence length for forward sequences is", ave_seq_length_f, "bp\n")

# same for reverse
if(paired_end){
  sampleRs <- if(length(fnRs)>=6) sample(fnRs,6) else fnRs[1:length(fnRs)]
  myRvsample <- ShortRead::readFastq(sampleRs)
  print(head(sread(myRvsample),10))
  print(tail(sread(myRvsample),10))
  ave_seq_length_r <- round(mean(sread(myRvsample)@ranges@width))
  cat("\naverage sequence length for reverse sequences is", ave_seq_length_f, "bp\n")
  ave_seq_length <- mean(c(ave_seq_length_f, ave_seq_length_r))
}

if(play_audio) beep(sound = sound_n)

rm(myFwsample)
if (paired_end) rm(myRvsample)
gc()
if(keep_time) toc()

#  check primers automatically --------------------------------------------

if(check_primers == "auto") {
  if(keep_time) tic("detecting primers...")
  cat("I will check for the occurrence of primers", "\n",
      "NOTE: this may take time, be patient...", "\n")
  
  # prefilter the sequences to remove ambiguous bases
  # check if a filtN directory is available, if not create
  cat("prefiltering sequences", "\n")
  if(!dir.exists(file.path("data","filtN"))) dir.create(file.path("data","filtN"))
  fnFs.filtN <- file.path("data", "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
  if(paired_end){
    fnRs.filtN <- file.path("data", "filtN", basename(fnRs))
  }
  if(paired_end){
    filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) 
  } else {
    filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = TRUE)
  }
  # count the occurrence of primers
  primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(ShortRead::readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
  }
  if(paired_end){
    primer_occ <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
                        FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
                        REV.ForwardReads = sapply(REV.orients, primerHits,fn = fnFs.filtN[[1]]), 
                        REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
  } else {
    primer_occ <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
                        REV.ForwardReads = sapply(REV.orients, primerHits,fn = fnFs.filtN[[1]])
    )
  }
  print(primer_occ)
  if(play_audio) beep(sound = sound_n)
  if(keep_time) toc()
  # sanity check; in the table  the forward primer should occur at the beginning of
  # forward sequences and its reverse complement should appear at the end of reverse sequences
  # reverse primer should appear at the beginning of reverse sequences and its RC at the
  # end of forward sequences
}  

#  use cutadapt -----------------------------------------------------------
cat("\n", "the use_cutadapt flag is set to ", use_cutadapt) 
cat("\n", "if you want to change the use_cutadapt do it below")
beep()
# no primers except for a few sequences
# use_cutadapt <- !use_cutadapt

# make sure you set the location of cutadapt correctly
if(use_cutadapt){
  if(keep_time) tic("primer removal")
  cat("removing primers...", "\n")
  # you can check if you have a cutadapt environment using reticulate package
  conda_list()
  # directory and file paths for processed files
  path.cut <- file.path("data", "cutadapt_out")
  if(!dir.exists(path.cut)) dir.create(path.cut)
  # get and save the project wd
  wd <- getwd()
  wd
  # need to do this because I will change the wd and then restore it
  if(exists("fnFs.cut")) rm(fnFs.cut) # needed to prevent problems if you rerun the script
  fnFs.cut <- file.path(wd, path.cut, basename(fnFs))
  if(paired_end) {
    if(exists("fnRs.cut")) rm(fnRs.cut)
    fnRs.cut <- file.path(wd, path.cut, basename(fnRs))
  }
  
  fnFs.filtN <- file.path(wd, fnFs.filtN)
  if(paired_end) fnRs.filtN <- file.path(wd, fnRs.filtN)
  # you need to set this to the directory containing your miniconda environments
  setwd(file.path("~","miniconda3"))
  getwd()
  cutadapt_loc <- file.path("envs", "cutadaptenv", "bin", "cutadapt")
  # if everything works this should return the version of cutadapt
  system2(cutadapt_loc, args = "--version")
  # create the reverse complements for the primers
  FWD.RC <- dada2:::rc(FWD)
  REV.RC <- dada2:::rc(REV)
  # Trim FWD and the reverse-complement of REV off of R1 (forward reads)
  R1.flags <- paste("-g", FWD, "-a", REV.RC) 
  # Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
  if(paired_end) R2.flags <- paste("-G", REV, "-A", FWD.RC) 
  # Run cutadapt
  for(i in seq_along(fnFs)) {
    if(paired_end){
      system2(cutadapt_loc, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                                     "-m", 20, # useful to remove 0 length sequences which cause problems afterwards
                                     "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                                     fnFs.filtN[i], fnRs.filtN[i])) # input files 
    } else {
      system2(cutadapt_loc, args = c(R1.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                                     "-m", 20, # useful to remove 0 length sequences which cause problems afterwards
                                     "-o", fnFs.cut[i],  # output files
                                     fnFs.filtN[i])) # input files   
    }
  }
  # check the occurrence of primers on the first sample
  if(paired_end){
    primer_occ_2 <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
                          FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
                          REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
                          REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
  } else {
    primer_occ_2 <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
                          REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]))
  }
  
  print(primer_occ_2)
  # Sanity check: a table with 0s should be returned
  
  # now reset working directory and 
  setwd(wd)
  # Forward and reverse fastq filenames have the format:
  cutFs <- sort(list.files(path.cut, pattern = "_1", full.names = TRUE))
  # in some very weird cases the file naming may be inconsistent with orientation
  if(!paired_end & !length(cutFs)==length(fnFs.cut)) {
    cutFs <- sort(list.files(path.cut, full.names = TRUE))
  }
  if(paired_end) cutRs <- sort(list.files(path.cut, pattern = "_2", full.names = TRUE))
  sampleFs <- if(length(cutFs)>=6) sample(cutFs,6) else cutFs[1:length(cutFs)]
  if(paired_end){
    sampleRs <- if(length(cutRs)>=6) sample(cutRs,6) else fnRs[1:length(cutRs)]
  }
  if(play_audio) beep(sound = sound_n)
  if(keep_time) toc()
}

# If you have PacBio data you should consider using dada2::removePrimers()
# if you want to see a report browseURL(report(qa(sampleFs))) (or replace with fnFs and fnRs)
# it may be useful if there is a mismatch in reads between forwards and reverse
# but it is slow

# creates quality profile plots --------------------------------

if(keep_time) tic("\ncreate quality profile plots")

# adapt this if you want to pick specific runs
# you should be able to determine from quality profiles if the quality scores are binned

plot_x_limit <- round(ave_seq_length/50)*50+25
toplot_fwd <- sampleFs
qplotfwd <- plotQualityProfile(toplot_fwd) +
  scale_x_continuous(limits = c(0,plot_x_limit), 
                     breaks = seq(0,plot_x_limit,25)) + 
  ggtitle("Fwd") + 
  theme(panel.grid.major.y = element_line(colour="grey75", linewidth = 0.5, linetype = 3),
        axis.text.x = element_text(angle = 90, hjust = 1))

qplotfwd 
# gray scale is the number of sequences of a given quality at a given position
# green average quality, continuous orange line median quality, 
# dashed 25th e 75th percentile;
# if the plot varies in length a continuous red line shows the % of sequences
# extending to this length

if(paired_end){
  toplot_rev <- sampleRs
  qplotrev <- plotQualityProfile(toplot_rev) + 
    scale_x_continuous(limits = c(0,plot_x_limit), 
                       breaks = seq(0,plot_x_limit,25)) + 
    ggtitle("Rev") + 
    theme(panel.grid.major.y = element_line(colour="grey75", linewidth=0.5, linetype = 3),
          axis.text.x = element_text(angle = 90, hjust = 1))
  qplotrev
}


# HINT make a note here of what you have done for reproducibility reasons: 
# primers have been removed

# save and remove object which won't be needed
if(paired_end){
  save(qplotfwd, qplotrev, file = str_c(Study,"qplots.Rdata"))
  rm(qplotfwd, qplotrev)
} else {
  save(qplotfwd, file = str_c(Study,"qplots.Rdata"))
  rm(qplotfwd)
}


# saving the workspace, with basic information
save.image(file = str_c(Study,".Rdata"))

if(keep_time) toc()
if(play_audio) beep(sound = sound_n)

# may be they have been quality processed

# filtering and trimming --------------------------------------------------

# creates directory data/filtered if it does not exist
if(!file_test("-d", filt_path)) dir.create(filt_path)

# creates file paths for filtered sequences

filtFs <- file.path(filt_path, basename(fnFs))

if(paired_end) {filtRs <- file.path(filt_path, basename(fnRs))}

# truncf position at which forward sequences will be truncated
# truncr idem for reverse
# quality score 30 1 error in 100, 40 1 in 10000, 20 1 in 100
# Novogene "clean" sequences are merged and quality processed, with primers removed
# reverse not applicable here: in some cases it is better to remove a few nt from the end
# NOTE: with ITS for fungi you should aim at a sum of 500 bp between forward and reverse to have good merging
truncf<- 0 # 0 in ITS
truncr<- 0 # 0 if not paired end and ITS
# the ITS DADA2 pipeline does not enforce a fixed length, but this might also result in removing too many
# sequences due to poor quality at the 3' end
trim_left = c(0,0) # use a length 2 vector c(x,y) if paired end or a single number if not
if(platform == "Ion_Torrent") trim_left <- trim_left+15
maxEEf = 2 # with very high quality data can be reduced to 1
maxEEr = 5 # 2 very restrictive, 5 does well in most cases, not needed for single end
trunc_q = 2 # with very high quality data can be increased up to 10-11
min_length <- 50
max_length <- 999 # not needed for Illumina and Ion Torrent, modify to max. exp. seq. length for F454
filter_and_trim_par <- as.data.frame(cbind(truncf, truncr, trim_left, 
                                           maxEEf, maxEEr, trunc_q, 
                                           max_length, min_length))

# matchIDs = true if prefiltered in QIIME;
if(use_cutadapt) {
  filt_sel <- which(basename(cutFs) %in% basename(fnFs))
  tofiltFs <- cutFs[filt_sel]
} else {
  tofiltFs <- fnFs
}
if(paired_end) {
  if(use_cutadapt) {
    filt_sel <- which(basename(cutRs) %in% basename(fnRs))
    tofiltRs <- cutRs[filt_sel]
  } else {
    tofiltRs <- fnRs
  }
}
# paired end 
out <- if(paired_end) {
  filterAndTrim(tofiltFs, filtFs, rev = tofiltRs, filt.rev = filtRs, 
                truncQ=trunc_q,
                trimLeft = trim_left,
                maxN=0, maxEE=c(maxEEf,maxEEr), 
                rm.phix=TRUE,  
                compress=TRUE, 
                minLen = min_length,
                multithread=TRUE) 
  # On Windows set multithread=FALSE
} else {
  # not paired end
  out <- filterAndTrim(tofiltFs, filtFs, rev = NULL, filt.rev = NULL, 
                       truncQ=trunc_q,
                       trimLeft = trim_left,
                       maxLen = max_length,
                       minLen = min_length,
                       maxN=0, maxEE=maxEEf, 
                       rm.phix=TRUE,  
                       compress=TRUE, 
                       multithread=TRUE) 
  # On Windows set multithread=FALSE
}
# a more complex alternative to optimize filter and trim parameters is to use
# Figaro https://github.com/Zymo-Research/figaro#figaro

head(out)
tail(out)
# summary info (use this to judge sequence loss and adjust your parameters)
minseq_left <- min(out[,2], na.rm = T)
maxseq_left <- max(out[,2], na.rm = T)
medseq_left <- median(out[,2], na.rm = T)
frac_loss <- (out[,1]-out[,2])/out[,1]
minfracloss <- min(frac_loss, na.rm = TRUE)
maxfracloss <- max(frac_loss, na.rm = TRUE)
medfracloss <- median(frac_loss, na.rm = TRUE)

if(play_audio) beep(sound = sound_n)

write_tsv(filter_and_trim_par, "filtertrimpars_conc.txt")

cat("After filtering there are between ", 
    minseq_left, 
    " and ", 
    maxseq_left, 
    " (median ", 
    medseq_left, 
    ") sequences left. The fraction of sequences lost after filtering is between ", 
    round(minfracloss, 2), 
    " and ", 
    round(maxfracloss, 2), 
    " (median ", 
    round(medfracloss, 2),
    ")\n",
    sep=""
    )

# better save workspace here
save.image(file = str_c(Study,".Rdata"))

# learning error rates ----------------------------------------------------

if(keep_time) tic("\nlearning error rates")

# as a default uses the first 1M sequences; if you want more nbases = xx
# another option could be randomize = T, to randomly pick samples for error
# estimates
# this may take VERY LONG, especially for quality binned data

qual_binned <- if_else(
  str_detect(platform, "miseq") | str_detect(platform, "Ion_Torrent") | str_detect(platform, "F454"),
  F, T)

# or else set manually
# qual_binned <- F

if(qual_binned){
  errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE, randomize = F)
  if(paired_end) errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE, randomize = F)
} else {
  errF <- learnErrors(filtFs, multithread=TRUE, randomize = F)
  if(paired_end) errR <- learnErrors(filtRs, multithread=TRUE, randomize = F)
  }

# plot error rates
plotErrors(errF, nominalQ=TRUE) + ggtitle("Fwd")
if(paired_end) plotErrors(errR, nominalQ=TRUE) + ggtitle("Rev")
# look at how the black line follows the points, should be monotonic
if(qual_binned){
  save(errF, file = str_c(Study,"_errF.Rdata"))
  new_errF_out <- getErrors(errF) %>%
    data.frame() %>%
    mutate_all(~case_when(. < X40 ~ X40,
                              . >= X40 ~ .)) %>% as.matrix()
  rownames(new_errF_out) <- rownames(getErrors(errF))
  colnames(new_errF_out) <- colnames(getErrors(errF))
  errF$err_out <- new_errF_out
  print(plotErrors(errF, nominalQ=TRUE) + ggtitle("Fwd_corr"))
  rm(new_errF_out)
  if(paired_end){
    save(errR, file = str_c(Study,"_errR.Rdata"))
    new_errR_out <- getErrors(errR) %>%
      data.frame() %>%
      mutate_all(~case_when(. < X40 ~ X40,
                            . >= X40 ~ .)) %>% as.matrix()
    rownames(new_errR_out) <- rownames(getErrors(errR))
    colnames(new_errR_out) <- colnames(getErrors(errR))
    errR$err_out <- new_errR_out
    print(plotErrors(errR, nominalQ=TRUE) + ggtitle("Rev_corr"))
    rm(new_errR_out)
  }
}

if(play_audio) beep(sound = sound_n)
if(keep_time) toc()

# dereplicate, infer ASV -------------------------------------------------------------

if(keep_time) tic("\ndereplicating, inferring ASVs, merging")

derepFs <- derepFastq(filtFs, verbose=TRUE)

if(paired_end) derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
if(paired_end) names(derepRs) <- sample.names


# inference of sequences; more accurate with pool = T
# acceptable run time with <1M sequences; becomes intractable with >10M
# an alternative is to use pool = "pseudo"

dadaFs <- if(platform == "Illumina") {
  dada2::dada(derepFs, err=errF, pool = F, multithread=TRUE)
} else {
  dada2::dada(derepFs, err=errF, pool = F, HOMOPOLYMER_GAP_PENALTY=-1, 
              BAND_SIZE=32, multithread=TRUE)
}

if(paired_end) dadaRs <- dada2::dada(derepRs, err=errR, pool = F, multithread=TRUE)

# better save workspace here
save.image(file = str_c(Study,".Rdata"))

if(keep_time) toc()
if(play_audio) beep(sound = sound_n)

# merge sequences ---------------------------------------------------------

# If paired end, when possible merge together the inferred forward and reverse sequences.
# with a good overlap (30 bp)
# set verbose = F if you want less output

# merge options
overlapping <- T # if F concatenation will be always performed
# use overlapping <- F if you loose all or most sequences after merging
minO = 25 # minimum overlap should be 20-30
maxM = 0 # maximum mismatch should be 0
# merge_option
# "fwd" will only use forward sequences (it is actually the same as operating on 
# non-paired end sequences
# "merge" will attempt a merge on all sequences
# "mixed" will try to merge all sequence it can and concatenate the others
# with this option no tree is returned
# minO and maxM are the minimum overlap and max mismatch
merge_option <- "mixed"

if(!paired_end | merge_option == "fwd"){
  mergers <- dadaFs
} else {
  if(keep_time) tic("\nmerging or concatenating paired ends")
  if (!overlapping){
    # when  fwd and rev are not able to cover the full amplicon length
    # with a given overlap justConcatenate = T adds 10 N between fwd and rev
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,
                          justConcatenate = T,
                          verbose=TRUE)
  } else {
    if(merge_option == "merge"){
      # does merging, can result in loss of longer sequences for ITS
      mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, 
                            maxMismatch = maxM,
                            minOverlap = minO,
                            verbose = TRUE)
    
      } else {
      mergers <- mixed_merge()
      }
    if(keep_time) toc()
  }
}

if(play_audio) beep(sound = sound_n)

# a bit of cleanup

# save and remove data not needed in further steps
if(paired_end){
  save(dadaFs, dadaRs, derepFs, derepRs, errF, errR,
       file = str_c(Study,"_dadaderep.Rdata"))
  rm(dadaRs, derepFs, derepRs, errF, errR)
} else {
  save(dadaFs, derepFs, errF,
       file = str_c(Study,"_dadaderep.Rdata"))
  rm(derepFs, errF)
}
save.image(file = str_c(Study,".Rdata"))
gc()

# build sequence table -----------------------------------------

if(keep_time) tic("\nmake sequence table")

seqtab.all <- makeSequenceTable(mergers)

# Inspect distribution of sequence lengths

lengthdistr <- table(nchar(getSequences(seqtab.all)))
# I should automate this step 
lengthdistr
barplot(lengthdistr, main = "sequence length distribution")

if(keep_time) toc()


# filter by size ----------------------------------------------------------

# ITS sequences are naturally variable in length and should not be filtered
# by length

# seqtab_f <- seqtab.all[,nchar(colnames(seqtab.all)) %in% seq(284, 472)]
seqtab_f <- seqtab.all
filt_seqs <- sum(seqtab_f)/sum(seqtab.all) 

# look at the abundance distribution
data.frame(nseqs = colSums(seqtab_f)) |>
  ggplot() +
  geom_histogram(mapping = aes(x= log10(nseqs))) +
  labs(
    title = "sequence length distribution", 
    x = "log10(seq_count)", 
    y = "number of individual sequences"
  )

# singletons will probably go away during chimera removal

save.image(file = str_c(Study,".Rdata"))

# remove bimera -------------------------------------------------------------
if(keep_time) tic("\nremove bimeras")
# may be very slow for large studies, you may want to check with a shorter version
# seqtab.nochim <- removeBimeraDenovo(seqtab_f[,1:10000])
cat("\nThe number of sequences prior to bimera removal is:", dim(seqtab_f)[2],"\n")

seqtab.nochim <- removeBimeraDenovo(seqtab_f, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
beep(sound=sound_n)
# better save workspace here
save.image(file = str_c(Study,".Rdata"))
# recheck length distribution and see if singletons are still there

data.frame(nseqs = colSums(seqtab.nochim)) |>
  ggplot() +
  geom_histogram(mapping = aes(x= log10(nseqs))) +
  labs(
    title = "sequence length distribution, post bimera removal", 
    x = "log10(seq_count)", 
    y = "number of individual sequences"
  )

# which ASVs are singletons or doubletons?
single_double <- which(colSums(seqtab.nochim)<=2)
length(single_double)/ncol(seqtab.nochim)
# fraction of singletons+doubletons here 2.5% (remember, these are ASVs, not OTUs) 
singletons <- which(colSums(seqtab.nochim)<=1)


doubletons <- which(colSums(seqtab.nochim)==2)
length(doubletons)

# If you want to remove singletons and doubletons set remove_s_d to T
remove_s_d <- F
seqtab.nochim.all <- seqtab.nochim
if(remove_s_d){
  seqtab.nochim <- seqtab.nochim[,-single_double]
}
cat("The number of sequences after bimera removal is:", dim(seqtab.nochim)[2], "\n")

# check the distribution of removals
summary(rowSums(seqtab.nochim)/rowSums(seqtab_f))

# sanity check ------------------------------------------------------------

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab_f), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

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

# View(track2)

track_long <- select(track2, sample_name, prop_filt:prop_nonchim) %>% 
  gather(key = "stage", value = "proportion", prop_filt:prop_nonchim) %>%
  mutate(stage = factor(stage, levels = c("prop_filt", "prop_denoised",
                                          "prop_merged", "prop_tabled",
                                          "prop_nonchim")))
track_long %>% ggplot(aes(x = sample_name, y = proportion, colour = stage)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0.5))
# very high losses for some samples
ggsave(filename = str_c("seqloss_",Study,".jpg"), width = 7, height = 5, 
       dpi = 150)
write_tsv(as.data.frame(track2), str_c("track_",Study,".txt"))

# further clean-up
save.image(file = str_c(Study,".Rdata"))
rm(dadaFs, mergers)
save.image(file = str_c(Study,"_small.Rdata"))

# assign taxonomy ---------------------------------------------------------

if(keep_time) tic("\nassign taxonomy")

taxdb_dir <- file.path("..","tax_db") # change this if the tax databases are elsewhere
list.files(taxdb_dir)

# assignment with UNITE
RC <- T # true by default

# UNITE
ref_fasta <- file.path(taxdb_dir, "sh_general_release_dynamic_s_04.04.2024.fasta")

# the most recent release is sh_general_release_dynamic_s_04.04.2024.fasta 

taxtab <- assignTaxonomy(seqtab.nochim, refFasta = ref_fasta, multithread = TRUE, 
                         tryRC = RC)

if(keep_time) toc()

if(play_audio) beep(sound=6)

# create a column with row names 

taxtab2 <- as_tibble(taxtab, rownames ="ASV")
if(!"Species" %in% colnames(taxtab)) taxtab2$Species <- NA_character_

# may be in the future add a variable with the Study, just to merge the database of ASVs
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
                                              ifelse(Species == "", Genus, paste(Genus, Species, sep =" "))))))))
  })


# as an alternative you can use DECIPHER::IdTaxa, which has been reported to
# be faster and more accurate than the naïve Bayesian classifier
# NEED TO CHECK IF THIS WORKS and provide the right database
# a number of training sets available here http://www2.decipher.codes/Downloads.html

useIdTaxa <- F
if(useIdTaxa){
  dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
  load(paste(taxdb_dir, "SILVA_SSU_r132_March2018.RData", sep="")) # CHANGE TO THE PATH OF YOUR TRAINING SET
  ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
  ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
  # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
  taxid <- t(sapply(ids, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
  }))
  colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
  taxtab2_IdTaxa <- taxid
}


# save workspace here: it may be truly large, you only need some of it
# remove unneeded objects
rm(ref_fasta)
save.image(file = str_c(Study,".Rdata"))

# save a smaller version of workspace which will be needed later
if(use_cutadapt){
  save(taxtab, seqtab.nochim, track2, Study, target, region, seq_accn, DOI,
       primer_f, primer_r, target1, target2, taxtab2, primer_occ, primer_occ_2,
       file = str_c(Study,"_small.Rdata"))
} else {
  save(taxtab, seqtab.nochim, track2, Study, target, region, seq_accn, DOI,
       primer_f, primer_r, target1, target2, taxtab2,
       file = str_c(Study,"_small.Rdata"))
}
if (play_audio) beep(sound = sound_n)


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


# Build phylogenetic tree ---------------------------------------------

# decipher can align only ~46k seqs
# https://www.bioconductor.org/packages/3.7/bioc/vignettes/DECIPHER/inst/doc/ArtOfAlignmentInR.pdf

# if this is >46k will crash
# but it is already prohibitive with 20k
# and usually not worth the effort

dotree <- F # avoid if overlapping == F

if(merge_option == "mixed" & dotree == T){
  cat(red("\nYour option for merging is 'mixed': do you really want to infer the phylogenetic tree?\n"))
}

if(dim(seqtab.nochim)[2]>10000 | overlapping == F) {
  dotree <- F
  }

if (dotree) {
  if(keep_time) tic("\nbuilding the phylogenetic tree")
  seqs <-
    dada2::getSequences(seqtab.nochim) # the collapse option is very interesting
  names(seqs) <- seqs # This propagates to the tip labels of the tree
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
  cat("creating distance matrix...","\n")
  dm <- phangorn::dist.ml(phang.align)
  # perform Neighbor joining
  cat("creating tree...","\n")
  treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order
  # internal maximum likelihood for tree
  #cat("estimating internal ML for tree...", "\n)
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

# check if any cols in seqtab have 0 sums
seq_sums <- rowSums(seqtab.nochim)
runs_to_keep <- which(seq_sums>0)
runs <- rownames(seqtab.nochim)[runs_to_keep]
# you may need to adapt this code
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
# you may need to adapt this code
if(data_type == "sra") {
  if(all(rownames(seqtab.nochim) %in% samdf$Run)){
    cat("\nsamples in sequence table match samples in metadata\n")
  } else {
    cat("\nWARNING samples in sequence table DO NOT match samples in metadata\n")
  }
} else {
  if(all(rownames(seqtab.nochim) %in% samdf$Library_Name)){
    cat("\nsamples in sequence table match samples in metadata\n")
  } else {
    cat("\nWARNING samples in sequence table DO NOT match samples in metadata\n")
  }
}
# rownames(seqtab.nochim) <- sapply(strsplit(rownames(seqtab.nochim), "_", fixed = T), `[`, 1)

# add final number of sequences and number of issues from track2
# you may need to adapt this code
if(data_type == "sra") {
  samdf <- left_join(samdf, select(track2, Run = sample_name, n_reads2 = nonchim, 
                                   n_issues))
} else {
  samdf <- left_join(samdf, select(track2, Library_Name = sample_name, n_reads2 = nonchim, 
                                   n_issues))
}

samdf <- as.data.frame(samdf)
if(data_type == "sra"){
  rownames(samdf) <- samdf$Run
} else {
  rownames(samdf) <- samdf$Library_Name
  }
# need ad hoc fix
samdf_copy <- samdf |>
  select(-description...23)
colnames(samdf_copy)[22]<-"description"
samdf_old <- samdf
samdf<-samdf_copy

# a transposed sequence table
seqtab_t = t(seqtab.nochim)

# combine in a phyloseq object
if(exists("fitGTR", mode = "list")){
  myphseq <- phyloseq(tax_table(taxtab), sample_data(samdf),
                      otu_table(seqtab_t, taxa_are_rows = TRUE),
                      phy_tree(fitGTR$tree))
} else {
  myphseq <- phyloseq(tax_table(taxtab), sample_data(samdf),
                      otu_table(seqtab_t, taxa_are_rows = TRUE))}


# save the object for further analysis
saveRDS(myphseq,str_c("data/", Study, "_ps.rds"))

# save the workspace
save.image(file = str_c(Study,"_small.Rdata"))

# prepare files for FMBN ----------------------------------------------

# The following commands are needed to prepare objects which, after further
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
# you need to set instrument and sequencing center manually if these are not in metadata
if(data_type == "sra"){
  instrument <- unique(samples$Instrument)[1]  # in some cases length >1
  seq_center <- unique(samples$Center_Name) 
  if(is_null(seq_center)) seq_center <- unique(samples$Center.Name)
} else {
  instrument <- "Illumina NovaSeq"
  seq_center <- "Novogene UK Ltd"
}
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

# information of geoloc (and names of the field) is very inconsistent:
# check the info in your sample metadata and adatp these commands
# use these if part or all of the geolocation information is missing
# samples$geo_loc_name_country <- "your country here"
# samples$geo_loc_name_country_continent <- "your continent here"
# samples$lat_lon <- NA_character_

# create label2 (to avoid numbers as first char.; s. can be removed later with
# tidyr::separate)

# need to be adjusted ad hoc
samples_copy <- samples
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
# on several days, yout need to reinitialize the log
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
                  primerr = primer_r,
                  primerfseq = FWD.orients,
                  primerrseq = REV.orients)
  log_print(primers, console = F)
  if(use_cutadapt){ 
    handling_primers <- list(
      checkprimers = check_primers,
      usecutatpt = use_cutadapt,
      primer_occ_pre = primer_occ,
      primer_occ_post = primer_occ_2
    )
  } else{
    handling_primers <- list(
      checkprimers = check_primers,
      usecutatpt = use_cutadapt
    )
  }
  log_print(handling_primers, console = F)
  log_print(filter_and_trim_par, console = F)
  if(paired_end) {merge_opt <- list(
    pend = paired_end,
    ovl = overlapping,
    mrgopt = merge_option,
    max_mismatch = maxM,
    min_ovlap = minO)} else {
      merge_opt <- list(
        pend = paired_end,
        ovl = overlapping,
        mrgopt = merge_option)
    }
  
  dotree
  log_close()
}


# Package citations -------------------------------------------------------

map(c(.cran_packages, .bioc_packages), citation)


# Credits and copyright ---------------------------------------------------

# Most of the script is taken from https://benjjneb.github.io/dada2/tutorial.html
# and https://benjjneb.github.io/dada2/ITS_workflow.html
# with several changes and adaptations

# Assume that this is overall under MIT licence

# Copyright 2024 Eugenio Parente
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
