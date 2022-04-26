################################################################################
# DADA2/Bioconductor pipeline for big data, modified
# part 1, split the fastq files in groups
#
# getHeader: a script to extract fastQ header v1.2 26/04/2022
################################################################################

# this script is designed to work in conjunction with a modified version of the
# DADA2 pipeline for big data (I am using it to manage studies with more than 
# 50 samples, but this would depend on the system you have)

# to use this script
# 1. copy the script to a new folder and create a RStudio project
# 2. create a new folder called data
# 3. inside the folder data create three folders: fastq, filtered, metadata
# 4. download from SRA the accession list and the study metadata for the
# any study you want to process (must be 16S) (here I am working on SRP59291)
# 5. using the sratoolkit download the fastq files from SRA and put them in the fastq folder

.cran_packages <- c("plyr", "tidyverse", "microseq", "stringr", "reshape2", "beepr")

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
# Load packages into session, and print package version
sapply(c(.cran_packages), require, character.only = TRUE)

set.seed(100)

sessionInfo()

platform <- "Illumina" # (or set to "Illumina" or "Ion_Torrent" or "F454")
paired_end <- T # set to true for paired end, false for single end

fastq_path <- file.path("data", "fastq")
filt_path <- file.path("data", "filtered")
# if(!file_test("-d", fastq_path)) {
#  dir.create(fastq_path)
# only needed if you want to create the directory and move the data from somewhere else


fns <- sort(list.files(fastq_path, full.names = TRUE))
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

# the size of the sample of sequences you want to use (see below)
n <-10
# to get a sample
sampleFs <- if(length(fnFs)>=n) sample(fnFs,n) else fnFs[1:length(fnFs)]
# to get all use fnFs

# a function to extract headers from fastq files
getHeader <- function(filename){
  myheader <- unlist(microseq::readFastq(filename)$Header)[1]
  dum<- cbind(header = myheader, run_name = basename(filename))
  return(dum)
}

# if you want to get sequences only on a sample replace fnFs with sampleFs
myFwsample <- ldply(fnFs, getHeader, .progress = "text")
beep(sound=6)

# needs to be optimized depending on the separators

identifiers <- as_tibble(myFwsample) %>%
  separate(header, into = c("sample", "machine_lane", "length"), 
           sep = " ", remove = F) 
# skip if necessary
skip_flag <- T
if(!skip_flag){
identifiers <- identifiers %>% separate(machine_lane, into = c("mach_lane_1", "mach_lane_2"),
                                        sep = "-", remove = F) %>%
  mutate(mach_lane_1 = as.factor(mach_lane_1))
  table(identifiers$mach_lane_1)
} else {
  identifiers <- identifiers %>% mutate(machine_lane = as.factor(machine_lane))
  table(identifiers$machine_lane)
}


# adapt this as needed, you can comfortably analyze 60-100 samples in one run
# make sure taht each group contains samples from a singli machine lane
identifiers$analysis_order <- c(rep(1,60), rep(2,55))




write_tsv(identifiers, file = "identifiers.txt")

# Package citations -------------------------------------------------------

map(.cran_packages, citation)





# Credits and copyright ---------------------------------------------------


# Assume that this is overall under MIT licence

# Copyright 2021, 2022 Eugenio Parente
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

