################################################################################
# DADA2/Bioconductor pipeline for big data, modified
# part 2, create sequence tables
#
# bigdata_part1_v8_3_1024_ITS
################################################################################

# This script is designed to process large studies using the
# DADA2 pipeline https://benjjneb.github.io/dada2/tutorial.html with options
# for large studies https://benjjneb.github.io/dada2/bigdata.html
# with modifications for ITS analysis https://benjjneb.github.io/dada2/ITS_workflow.html
# a further script (bigdata_part2) will carry out further processing needed to ready the output for import
# into FoodMicrobionet

# This version of the script (v8.1, 7/8/24) includes options
# for single end/paired end data sets obtained with Illumina or 454 or
# Ion Torrent and will perform steps up to the creation of the sequence table
# including optional primer removal with cutadapt
# In addition it includes a workaround for Illumina Novaseq sequences described here:
# https://github.com/benjjneb/dada2/issues/791 and a workaround for  handling
# ITS sequences which do not merge properly due to excess length described here:
# https://github.com/benjjneb/dada2/issues/537
# In this version UNITE general release for fungi is used as a taxonomic reference
# For each iteration of the script over the identifiers.txt table, created with
# the getHeader script, a new sequence table is created and saved from the
# corresponding group fo sequences. These are then assembled and processed
# usign the bigdata_part2 script

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
#    4.1 put the folder you received from Novogene in the data folder
#    4.2 create the filtered and metadata folders within the data folder
#    4.3 put the metadata in the metadata folder; the first column must match 
#        sample names (as received from Novogene, these are also the names of 
#        the folders containing the sequences)
# 5. download the most recent taxonomy reference files from 
#    https://benjjneb.github.io/dada2/training.html
#    and put them in a folder called taxdb at the same level of the directory 
#    containing your project (will be needed in the second part of this script)
# 6. an optional table with commonly used primer pairs (primer_pairs_fungi.txt) 
#    can be downloaded from our repository and is used by a function to 
#    check primers; it should be placed in the folder containing the folder of 
#.   your project and the taxonomic reference folder.
# 7. make sure you have all the information you need (primers, platform, region, etc.)
# You can easily adapt the script to your own data
# Make sure you have run the getHeader script and created the identifiers.txt table
# Set options for group and metadata in the relevant sections
# If you want to use cutadapt to remove primers install it using miniconda
# A good tutorial is here https://astrobiomike.github.io/unix/conda-intro
# You can easily adapt the script to your own data
# Run the script one section at a time and check the output and warnings, if any

# load packages -----------------------------------------------------------

.cran_packages <- c("tidyverse", "parallel", "beepr", "tictoc", "reticulate", 
                    "logr", "R.utils")
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

set.seed(100L)

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

# standard options for handling primers: you have to change this manually 
# if you want to do otherwise
check_primers <- "auto" # "visually" if you want to do visually, "auto" otherwise
use_cutadapt <- T # use cutadapt to remove primers, takes time and space, use wisely
use_primer_table <- T # uses a primer table to check if the primer sequences are OK
# put the location of your primer table inside this if
if(use_primer_table){
  primer_table_path <- file.path("..", "primer_pairs_fungi.txt")
}


# functions ---------------------------------------------------------------

# this will load three functions which are going to be used in the pipeline
# this file should be in the folder containing your project directory
# otherwise you need to change the path
source(file.path("..","bioconductor_pip_ITS_functions.R"))

# set group ---------------------------------------------------------------

# open identifiers (a file containing the grouping of your runs)

mygroup <- 3
# needs to be adapted for each iteration
# remember there is also a copy created which is never modified
if (mygroup==1){
  identifiers <- read_tsv("identifiers_copy.txt")
  identifiers <- identifiers %>% 
    mutate(run_2 = str_replace(run_name, "_1", "_2")) %>%
    select(header, analysis_order, forward = run_name, reverse = run_2)
  write_tsv(identifiers, file = "identifiers.txt")
} else {
  identifiers <- read_tsv("identifiers.txt")
}

# opts_chunk$set(cache = FALSE,fig.path="dadafigure/")
# read_chunk(file.path("src", "bioinformatics.R"))
# creating information for the study and sample dataframes

if (mygroup==1){
  # change this as appropriate
  data_type <- "sra" 
  # alternatives are "sra", for data downloaded from sra
  # "novogene_raw", for data obtained from novogene, with primer not removed
  # "novogene_clean", for data obtained from novogene, data with primer removed
  
  # creating information for the study and sample data frames
  # when using data other than those downloaded from SRA replace
  # the accession number with whichever identifier you want to use
  
  Study <- "SRP235327" # PRJNA594469
  target <- "ITS region"
  region <- "ITS1"
  seq_accn <- Study
  DOI <- "10.1128/msphere.00534-20"
  
  # information on the platform and arrangement
  
  platform <- "Illumina_miseq" #  (set to "Illumina_miseq", "Illumina_novaseq", "Illumina_HiSeq" or "Ion_Torrent" or "F454")
  paired_end <- T # set to true for paired end, false for single end
  overlapping <- T # here set to F because of poor quality of reverse sequences which prevented merging
  
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
  primer_f <- "ITS1F" 
  primer_r <- "ITS2"
  
  # NOTE: be extra careful in indicating primer sequences because this will affect
  # primer detection and primer removal by cutadapt
  
  FWD <- "TTGGTCATTTAGAGGAAGTAA"  ## CHANGE THIS to your forward primer sequence
  REV <- "GCTGCGTTCTTCATCGATGC" ## CHANGE THIS to your reverse primer sequence
  
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
  
  FWD.orients <- allOrients(FWD)
  REV.orients <- allOrients(REV)
  print(FWD.orients)
  print(REV.orients)
  
  # save general data
  gen_data <- list(Study = Study, target = target, region = region, DOI = DOI,
                   data_type = data_type,
                   platform = platform, overlapping = overlapping, paired_end = paired_end,
                   primer_f = primer_f, primer_r = primer_r, 
                   primerf_seq = FWD, primerr_seq = REV,
                   target1 = target1)
  saveRDS(gen_data, file = "gen_data.RDS")
  } else {
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
  FWD.orients <- allOrients(FWD)
  REV.orients <- allOrients(REV)
  print(FWD.orients)
  print(REV.orients)
}

# create subdirectories ----------------------------------------------------------
# creates subdirectories and get names of the fastq files to process
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
      # if you want the seqs without primer and adapter
      # file_index <- which(!str_detect(file_list[[i]], "raw") & !str_detect(file_list[[i]], "Frags"))
      file_index <- which(str_detect(file_list[[i]], "raw"))
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

# extracts the selection for the specific group being processed
my_selection <- identifiers %>%
  filter(analysis_order == mygroup)
group_to_extract <- c(
  pull(my_selection, forward),
  pull(my_selection, reverse)
)
indices <- which(basename(fns) %in% group_to_extract)

fns <- fns[indices]

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
if(data_type == "novogene_raw"){
  sample.names <- str_remove(sample.names, "\\.raw")
}


# early check for sequences and metadata --------------------------------

# path to the metadata file (needs to be adapted)
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

# check for occurrence of primers and adapters on a sample of forward and
# reverse sequences, done only for the first group (if you suspect there are differences
# do it on multiple groups)
if (mygroup == 1) {
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
  if(keep_time) toc()
}
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

if(mygroup == 1){
  if(keep_time) tic("\ncreate quality profile plots")
  
  # adapt this if you want to pick specific runs
  
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
  
  rm(myFwsample)
  if (paired_end) rm(myRvsample)
  
  # HINT make a note here of what you have done for reproducibility reasons: 
  # primers have been removed
  
  # save and remove object which won't be needed
  if(paired_end){
    save(qplotfwd, qplotrev, file = str_c(gen_data$Study,"_", mygroup, "qplots.Rdata"))
    rm(qplotfwd, qplotrev)
  } else {
    save(qplotfwd, file = str_c(gen_data$Study,"_", mygroup, "qplots.Rdata"))
    rm(qplotfwd)
  }
}


# saving the workspace, with basic information
save.image(str_c(gen_data$Study,"_", mygroup, ".Rdata"))

if(keep_time) toc()
if(play_audio) beep(sound = sound_n)


# sequence filtering ------------------------------------------------------

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
truncf<- NULL # NULL in ITS
truncr<- NULL # NULL if not paired end and ITS
# the ITS DADA2 pipeline does not enforce a fixed length, but this might also result in removing too many
# sequences due to poor quality at the 3' end
trim_left = 0 # use a length 2 vector c(x,y) if paired end or a single number if not
if(platform == "Ion_Torrent") trim_left <- trim_left+15
maxEEf = 2 # with very high quality data can be reduced to 1
maxEEr = 2 # 2 very restrictive, 5 does well in most cases, not needed for single end
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


if (mygroup == 1) write_tsv(filter_and_trim_par, "filtertrimpars_conc.txt")
save.image(str_c(gen_data$Study,"_", mygroup, ".Rdata"))

# learn error rates ------------------------------------------

if(keep_time) tic("\nlearning error rates")

# as a default uses the first 1M sequences; if you want more nbases = xx
# another option could be randomize = T, to randomly pick samples for error
# estimates
# this may take VERY LONG, especially for novaseq

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
  save(errF, file = str_c(gen_data$Study,"_", mygroup,"_errF.Rdata"))
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
    save(errR, file = str_c(gen_data$Study,"_", mygroup, "_errR.Rdata"))
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

save.image(str_c(gen_data$Study,"_", mygroup, ".Rdata"))
if(play_audio) beep(sound = sound_n)
if(keep_time) toc()

# infer ASVs and merge  --------------------------------------------------------------

if(keep_time) tic("\ninferring AVSs and merging")

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

if(mygroup == 1){
  gen_data$merge_option <- merge_option
  saveRDS(gen_data, file = "gen_data.RDS")
}

# infer ASVs and merge if applicable
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
n_samples <- length(sample.names)
getN <- function(x) sum(getUniques(x)) # a function for getting number of sequences out of dada
denoised <- vector(mode = "integer", length = n_samples)

# ad hoc if you need to remove a problematic sample, here SRR10613791
# which(sample.names == "SRR10613791")
# length(sample.names)
# sample.names <- sample.names[-40]
# length(sample.names)

for(sam in sample.names) {
  sample_index_F <- which(str_detect(filtFs, sam))
  cat("Processing:", sam, "-", sample_index_F,
      "of", n_samples, "\n")
  derepF <- derepFastq(filtFs[[sample_index_F]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  denoised[sample_index_F] <- getN(ddF)
  if(paired_end){
    sample_index_R <- which(str_detect(filtRs, sam))
    derepR <- derepFastq(filtRs[[sample_index_R]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
  }
  if(!paired_end | merge_option == "fwd"){
    merger <- dadaFs
  } else {
    if(keep_time) tic("\nmerging or concatenating paired ends")
    if (!overlapping){
      # when  fwd and rev are not able to cover the full amplicon length
      # with a given overlap justConcatenate = T adds 10 N between fwd and rev
      merger <- mergePairs(ddF, derepF, ddR, derepR,
                            justConcatenate = T,
                            verbose=TRUE)
    } else {
      if(merge_option == "merge"){
        # does merging, can result in loss of longer sequences for ITS
        merger <- mergePairs(ddF, derepF, ddR, derepR, 
                              maxMismatch = maxM,
                              minOverlap = minO,
                              verbose = TRUE)
        
      } else {
        merger <- mixed_merge_bigdata()
      }
      if(keep_time) toc() 
    }
  }

  mergers[[sam]] <- merger
}
lapply(mergers, class)


if(play_audio) beep(sound = sound_n)
if(keep_time) toc()

rm(derepF)
if(paired_end) rm(derepR)
gc()

# Build sequence table  -----------------------------------------------
# remove if the object is NULL
null_mergers <-which(sapply(mergers, class)=="NULL")

seqtab <- makeSequenceTable(mergers[-null_mergers])

# save data
save.image(file = str_c(gen_data$Study,"_", mygroup, ".Rdata"))

# Inspect distribution of sequence lengths

lengthdistr <- table(nchar(getSequences(seqtab)))
# I should automate this step
lengthdistr
barplot(lengthdistr)
# to remove sequences much shorter or longer than expected

# filter by size ----------------------------------------------------------

# ITS sequences are naturally variable in length and should not be filtered
# by length

# seqtab_f <- seqtab.all[,nchar(colnames(seqtab.all)) %in% seq(284, 472)]
seqtab_f <- seqtab
filt_seqs <- sum(seqtab_f)/sum(seqtab) 

# look at the abundance distribution
data.frame(nseqs = colSums(seqtab_f)) |>
  ggplot() +
  geom_histogram(mapping = aes(x= log10(nseqs))) +
  labs(
    title = "sequence length distribution", 
    x = "log10(seq_count)", 
    y = "number of individual sequences"
  )

# rename and save sequence table
assign(str_c("seqtab_f_",mygroup), seqtab_f)
save(seqtab_f, file = str_c("seqtab_f_", mygroup, ".Rdata"))

# sanity check (without chimera step)

track <- cbind(out[-null_mergers,], denoised[-null_mergers], 
               sapply(mergers[-null_mergers], getN), 
               rowSums(seqtab_f))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled")
rownames(track) <- sample.names
head(track)

saveRDS(track, file = str_c("track_",mygroup,".RDS"))
save.image(file = str_c(Study,"_",mygroup,".Rdata"))


# optionally save a log ---------------------------------------------------

if(use_logr){
  log_path <- file.path(str_c(Study, mygroup, "log", sep = "_"))
  log_open(log_path)
  gen_options <- list(
    audio = play_audio,
    sound = sound_n,
    timing = keep_time,
    verbose = verbose_output,
    group = mygroup
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
  log_close()
}


# Package citations -------------------------------------------------------

map(c(.cran_packages, .bioc_packages), citation)


# Credits and copyright ---------------------------------------------------

# Most of the script is taken from https://benjjneb.github.io/dada2/tutorial.html
# or https://benjjneb.github.io/dada2/bigdata.html
# with some changes and adaptations, mostly to allow the automation of the
# second step (sequence merging, chimera removal, taxonomic assignment etc.)

# Assume that this is overall under MIT licence

# Copyright 2022, 2024 Eugenio Parente
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

