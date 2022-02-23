################################################################################
# DADA2/Bioconductor pipeline for big data, modified
# part 2, create sequence tables
#
# bigdata_part1_v5_1021
################################################################################

# This script is designed to process large studies using the
# DADA2 pipeline https://benjjneb.github.io/dada2/tutorial.html with options
# for large studies https://benjjneb.github.io/dada2/bigdata.html
# a further script (bigdata_part2) will carry out further processing needed to ready the output for import
# into FoodMicrobionet

# This version of the script (v5_0, 10/2021) includes options
# for single end/paired end data sets obtained with Illumina or 454 or
# Ion Torrent
# and will perform steps up to teh creation of the sequence table
# For each iteration of the script over the identifiers.txt table, created with
# the getHeader script, a new sequence table is created and saved from the
# corresponding group fo sequences. These are then assembled and processed
# usign the bigdata_part2 script

# to use this script
# 1. copy the script to a new folder and create a RStudio project
# 2. create a new folder called data
# 3. inside the folder data create three folders: fastq, filtered, metadata
# 4. download from SRA the accession list and the study metadata for the
# any study you want to process (must be 16S) (here I am working on SRP229651)
# 5. using the sratoolkit download the fastq files from SRA and put them in the fastq folder
# 6. put the metadata i the metadata folder
# 7. download the taxonomy reference files from https://benjjneb.github.io/dada2/training.html
# and put them in a folder called taxdb at the same level of the directory containing your project
# make sure you have all the information you need (primers, platform, region, etc.)
# You can easily adapt the script to your own data
# Make sure you have run the getHeader script and created the identifiers.txt table
# Set options for group (line 71) and metadata (lines 84-128)
# Run the script one section at a time

# load packages -----------------------------------------------------------

.cran_packages <- c("tidyverse", "gridExtra", "knitr", "stringr", "reshape2", "beepr")
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

set.seed(100)

sessionInfo()

# set group ---------------------------------------------------------------

# open identifiers (a file containing the grouping of your runs)
identifiers <- read_tsv("identifiers.txt")
mygroup <- 2
# needs to be adapted
if (mygroup==1){
  identifiers <- identifiers %>% 
    mutate(run_2 = str_replace(run_name, "_1", "_2")) %>%
    select(header, analysis_order, forward = run_name, reverse = run_2)
  write_tsv(identifiers, file = "identifiers.txt")
}

# opts_chunk$set(cache = FALSE,fig.path="dadafigure/")
# read_chunk(file.path("src", "bioinformatics.R"))
# creating information for the study and sample dataframes

if (mygroup==1){
  Study <- "SRP229651"
  target <- "16S RNA gene"
  region <- "V1-V3"
  seq_accn <- Study
  DOI <- "10.1111/1462-2920.15407"
  
  # information on the platform and arrangement
  
  platform <- "Illumina" # (or set to "Illumina" or "Ion_Torrent" or "F454")
  paired_end <- T # set to true for paired end, false for single end
  overlapping <- F # here set to F because of poor quality of reverse sequences which prevented merging
  
  
  # Frequently used primer sets: --------------------------------------------
 
  # V1-V3 
  # Gray28F (5′-TTTGATCNTGGCTCAG) 16 bp and Gray519r (5′-GTNTTACNGCGGCKGCTG) 18 bp, 509 bp
  # 27F 5′ AGAGTTTGATCMTGGCTCAG3’ — 519R 5′ GWATTACCGCGGCKGCTG3′
  # 28F GAGTTTGATCNTGGCTCAG 19 bp 
  
  # 
  # V3-V4
  # S-D-Bact-0341-b-S-17/S-d-Bact-0785-a-A-21 primer pair with an amplicon size of 464 bp (Klindworth et al., 2012; 
  # also known as Bakt_341F/Bakt_805R primer pair designed by Sinclair et al. (2015)
  # forward CCTACGGGNGGCWGCAG 17 bp 
  # reverse GACTACHVGGGTATCTAATCC 21 bp
  # Uni340F (5′ CCTACGGGRBGCASCAG-3′) 17 bp and Bac806R (5′GGACTACYVGGGTATCTAAT 3′) 20 bp
  # 357 F (5’-CTCCTACGGGAGGCAGCAG-3’) 19 bp and 939R (5’-CTTGTGCGGGCCCCCGTCAATTC-3’) 23 bp - 605 bp
  #
  # V4
  # 515F (5'- GTGCCAGCMGCCGCGGTAA-3' ) 19 bp and 806R (5' -GGACTACHVGGGTWTCTAAT-3’) 20 bp length  311 bp
  # 515 F (5′- GTGBCAGCMGCCGCGGTAA - 3′) (Hugerth et al., 2014) and Pro--mod-805 R (5′-GACTACNVGGGTMTCTAATCC - 3′) 21 bp length 310 bp
  #
  # V4-V5 (515F, 5′-GTGCCAGCMGCCGCGGTAA-3′ 19 bp; 926R, 5′-CCGTCAATTCMTTTRAGT-3′  18 bp 429 bp) 
  
  
  # expected amplicon length ? including primers
  primer_f <- "28F" 
  primer_r <- "Gray519r"  
  target1 <- "16S_RNA"
  target2 <- region
  
  # V3-V4
  # forward primer 17 bp
  # reverse primer 21 bp
  # expected length 464
  
  # save general data
  gen_data <- list(Study = Study, target = target, region = region, DOI = DOI,
                   platform = platform, overlapping = overlapping, paired_end = paired_end,
                   primer_f = primer_f, primer_r = primer_r, target1 = target1)
  saveRDS(gen_data, file = "gen_data.RDS")
} else {
  gen_data <- readRDS(file = "gen_data.RDS")
  Study <- gen_data$Study
  target <- gen_data$target
  region <- gen_data$region
  seq_accn <- Study
  DOI <- gen_data$DOI
  
  # information on the platform and arrangement
  
  platform <- gen_data$platform 
  overlapping <- gen_data$overlapping
  paired_end <- gen_data$paired_end 
  primer_f <- gen_data$primer_f 
  primer_r <- gen_data$primer_r 
  target1 <- gen_data$target1
  target2 <- region
}



# create subdirectory ----------------------------------------------------------
# subdirectory data must already be in the wd

fastq_path <- file.path("data", "fastq")
filt_path <- file.path("data", "filtered")
# if(!file_test("-d", fastq_path)) {
#  dir.create(fastq_path)
# only needed if you want to create the directory and move the data from somewhere else


# get all filenames
fns <- sort(list.files(fastq_path, full.names = TRUE))

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
sample.names <- case_when(
  all(str_detect(basename(fnFs),"_")) ~ sapply(strsplit(basename(fnFs), "_", fixed = T), `[`, 1),
  all(str_detect(basename(fnFs),".")) ~ sapply(strsplit(basename(fnFs), ".", fixed = T), `[`, 1),
)

# check for occurrence of primers and adaptors on a sample of forward and
# reverse sequences, done only for the first lane (if you suspect there are differences
# do it on multiple groups)
if (mygroup == 1) {
  sampleFs <- if(length(fnFs)>=6) {sample(fnFs,6)} else {fnFs[1:length(fnFs)]}
  myFwsample <- ShortRead::readFastq(sampleFs)
  print(head(sread(myFwsample),10))
  print(tail(sread(myFwsample),10))
  ave_seq_length <- round(mean(sread(myFwsample)@ranges@width))
  
  if(paired_end){
    sampleRs <- if(length(fnRs)>=6) sample(fnRs,6) else fnRs[1:length(fnRs)]
    myRvsample <- ShortRead::readFastq(sampleRs)
    print(head(sread(myRvsample),10))
    print(tail(sread(myRvsample),10))
  }
  
  # sequences contain primers at the 5' end 
  rm(myFwsample)
  if (paired_end) rm(myRvsample)
  
  # creates quality profile plots --------------------------------
  # adapt this if you want to pick specific runs
  plot_x_limit <- round(ave_seq_length/100)*100+25
  toplot_fwd <- sampleFs
  qplotfwd <- plotQualityProfile(toplot_fwd) + 
    scale_x_continuous(limits = c(0,plot_x_limit), 
                       breaks = seq(0,plot_x_limit,25)) + 
    ggtitle("Fwd") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  print(qplotfwd)
  # grayscale is the number of sequences of a given quality at a given position
  # green average quality
  # continuous orange line median quality, dashed 25th e 75th percentile
  # if the plot varies in length a continuous line shows the % of sequences
  # extending to this length
  
  if(paired_end){
    toplot_rev <- sampleRs
    qplotrev <- plotQualityProfile(toplot_rev) + 
      scale_x_continuous(limits = c(0,plot_x_limit), 
                         breaks = seq(0,plot_x_limit,25)) + 
      ggtitle("Rev") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(qplotrev)
  }

  
  # save and remove object which won't be needed
  if(paired_end){
    save(qplotfwd, qplotrev, file = str_c(Study,"qplots_",mygroup,"_.Rdata"))
    rm(qplotfwd, qplotrev)
  } else {
    save(qplotfwd, file = str_c(Study,"qplots_",mygroup,".Rdata"))
    rm(qplotfwd)
  }
  
  
  # saving the workspace
  save.image(file = str_c(Study,"_", mygroup, ".Rdata"))
  beep(sound=6)
}

# sequence filtering ------------------------------------------------------

# creates directory data/filtered if it does not exist
if(!file_test("-d", filt_path)) dir.create(filt_path)
# creates file paths for filtered sequences
filtFs <- file.path(filt_path, basename(fnFs))

if(paired_end) filtRs <- file.path(filt_path, basename(fnRs))


# truncf position at which forward sequences will be truncated
# truncr idem for reverse
# quality score 30 1 error in 100, 40 1 in 10000, 20 1 in 100

truncf<-290
truncr<-235 
trim_left = c(19,18) 
if(platform == "Ion_Torrent") trim_left <- trim_left+15
maxEEf = 2 # with very high quality data can be reduced to 1
maxEEr = 5 # 2 very restrictive, 5 does well in most cases, not needed for single end
trunc_q = 2 # with very high quality data can be increased up to 10-11
max_length <- 999 # not needed for Illumina and Ion Torrent, modify to max. exp. seq. length for F454
filter_and_trim_par <- as.data.frame(cbind(truncf, truncr, trim_left, 
                                           maxEEf, maxEEr, trunc_q, max_length))

# matchIDs = true if prefiltered in QIIME;
# paired end

out <- if(paired_end) {
  filterAndTrim(fnFs, filtFs, rev = fnRs, filt.rev = filtRs, 
                truncQ=trunc_q,
                truncLen=c(truncf,truncr),
                trimLeft = trim_left,
                maxN=0, maxEE=c(maxEEf,maxEEr), 
                rm.phix=TRUE,  
                compress=TRUE, 
                multithread=TRUE) 
  # On Windows set multithread=FALSE
} else {
  # not paired end
  max_length <- ifelse(max_length>truncf, Inf, max_length)
  out <- filterAndTrim(fnFs, filtFs, rev = NULL, filt.rev = NULL, 
                       truncQ=trunc_q,
                       truncLen=truncf,
                       trimLeft = trim_left,
                       maxLen = max_length,
                       maxN=0, maxEE=c(maxEEf), 
                       rm.phix=TRUE,  
                       compress=TRUE, 
                       multithread=TRUE) 
  # On Windows set multithread=FALSE
  }


head(out)
tail(out)
beep(sound=6)
# some loss but number of sequences acceptable
# summary info (use this to judge sequence loss and adjust your parameters)
minseq_left <- min(out[,2], na.rm = T)
maxseq_left <- max(out[,2], na.rm = T)
medseq_left <- median(out[,2], na.rm = T)
frac_loss <- (out[,1]-out[,2])/out[,1]
minfracloss <- min(frac_loss, na.rm = TRUE)
maxfracloss <- max(frac_loss, na.rm = TRUE)
medfracloss <- median(frac_loss, na.rm = TRUE)

cat("After filtering there are between", 
    minseq_left, 
    "and", 
    maxseq_left, 
    "(median", 
    medseq_left, 
    ") sequences left. The fraction of remaining sequences after filtering is between", 
    round(minfracloss, 2), 
    "and", 
    round(maxfracloss, 2), 
    "(median", 
    round(medfracloss, 2),
    ")\n"
)

if (mygroup == 1) write_tsv(filter_and_trim_par, "filtertrimpars_conc.txt")

# learn error rates, infer ASVs, make sequence table ---------------------------

errF <- learnErrors(filtFs, multithread=TRUE, randomize = F)
if(paired_end) errR <- learnErrors(filtRs, multithread=TRUE, randomize = F)

# plot error rates
plotErrors(errF, nominalQ=TRUE) + ggtitle("Fwd")
if(paired_end) plotErrors(errR, nominalQ=TRUE) + ggtitle("Rev")
# look at how the black line follows the points

# save data
save.image(file = str_c(gen_data$Study,"_", mygroup, ".Rdata"))
beep(sound=6)
# infer ASVs and merge if applicable
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
n_samples <- length(sample.names)

getN <- function(x) sum(getUniques(x)) # a function for getting number of sequences out of dada

denoised <- vector(mode = "integer", length = n_samples)
for(sam in sample.names) {
  sample_index_F <- which(str_detect(filtFs, sam))
  cat("Processing:", sam, "-", sample_index_F,
      "of", n_samples, "\n")
  derepF <- derepFastq(filtFs[[sample_index_F]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  denoised[sample_index_F] <- getN(ddF)
  sample_index_R <- which(str_detect(filtRs, sam))
  derepR <- derepFastq(filtRs[[sample_index_R]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  if(paired_end){
    minO = 20 # should be 20-30
    maxM = 0 # should be 0
    # justConcatenate = T adds 10 N between fwd and rev
    if (overlapping){
      merger <- mergePairs(ddF, derepF, ddR, derepR, 
                            maxMismatch = maxM,
                            minOverlap = minO,
                           verbose=TRUE)
    } else {
      merger <- mergePairs(ddF, derepF, ddR, derepR,
                            justConcatenate = T,
                           verbose=TRUE)
    }
  } else {
    merger <- ddF
  }
  mergers[[sam]] <- merger
}
beep(sound=6)

rm(derepF); rm(derepR)


# Construct sequence table 
seqtab <- makeSequenceTable(mergers)

# save data
save.image(file = str_c(gen_data$Study,"_", mygroup, ".Rdata"))

# Inspect distribution of sequence lengths

lengthdistr <- table(nchar(getSequences(seqtab)))
# I should automate this step
lengthdistr
barplot(lengthdistr)
# to remove sequences much shorter or longer than expected

# should be the same used in the first group
seqtab_f <- seqtab[,nchar(colnames(seqtab)) %in% seq(497,499)]
sum(seqtab_f)/sum(seqtab)
# 100% seqs remaining after filtering
# look at the abundance distribution
qplot(log10(colSums(seqtab_f)))
# several singletons, will probably go away during chimera removal

# rename and save sequence table
assign(str_c("seqtab_f_",mygroup), seqtab_f)
save(seqtab_f, file = str_c("seqtab_f_", mygroup, ".Rdata"))

# sanity check (without chimera step)

track <- cbind(out, denoised, 
               sapply(mergers, getN), 
               rowSums(seqtab_f))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled")
rownames(track) <- sample.names
head(track)

saveRDS(track, file = str_c("track_",mygroup,".RDS"))
save.image(file = str_c(Study,"_",mygroup,".Rdata"))

# Credits and copyright ---------------------------------------------------

# Most of the script is taken from https://benjjneb.github.io/dada2/tutorial.html
# or https://benjjneb.github.io/dada2/bigdata.html
# with some changes and adaptations, mostly to allow the automation of the
# second step (sequence merging, chimera removal, taxonomic assignment etc.)

# Assume that this is overall under MIT licence

# Copyright 2021 Eugenio Parente
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

