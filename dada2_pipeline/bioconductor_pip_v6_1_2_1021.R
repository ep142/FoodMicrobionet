################################################################################
# DADA2/Bioconductor pipeline, modified, v6_1_3, 22/11/2021
################################################################################

# This script is designed to process reasonably large studies using the
# DADA2 pipeline https://benjjneb.github.io/dada2/tutorial.html
# and to carry out further processing needed to ready the output for import
# into FoodMicrobionet

# This version of the script (v6_0, 09/2021) includes options
# for single end/paired end data sets obtained with Illumina or 454 or
# Ion Torrent
# In addition Greengenes is not used anymore for taxonomic assignment
# as it is no longer being maintained, and SILVA v138.1 replaces v132 and v138
# Finally, the script assembles an object which can be used in the future for
# redoing the taxonomic assignment in an automated way

# to use this script
# 1. copy the script to a new folder and create a RStudio project
# 2. create a new folder called data
# 3. inside the folder data create three folders: fastq, filtered, metadata
# 4. download from SRA the accession list and the study metadata for the
# any study you want to process (must be 16S) 
# 5. using the sratoolkit download the fastq files from SRA and put them in the fastq folder
# 6. put the metadata in the metadata folder
# 7. download the taxonomy reference files from https://benjjneb.github.io/dada2/training.html
# and put them in a folder called taxdb at the same level of the directory containing your project
# make suere you have all the information you need (primers, platform, region, etc.)
# You can easily adapt the script to your own data
# Run the script one section at a time

# Install/load packages ---------------------------------------------------

.cran_packages <- c("tidyverse", "parallel","gridExtra", "knitr", "stringr", 
                    "phylotools", "beepr", "tictoc")
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
R.Version() # in the future may me check that the version running is compatible
# with the script (it only matters if somebody other than myself is using the script)
sessionInfo()

# opts_chunk$set(cache = FALSE,fig.path="dadafigure/")
# read_chunk(file.path("src", "bioinformatics.R"))

# other setup operations
opar <- par(no.readonly=TRUE) 
par(ask=F) 
# set.seed(1234) 
# play audio notifications
play_audio <- T
# time important steps
keep_time <- T
# verbose output: will print additional objects and messages
verbose_output <- T

if(play_audio) beep(sound = 6) # notification

# creating information for the study and sample data frames

Study <- "SRP259291"
target <- "16S RNA gene"
region <- "V3-V4"
seq_accn <- Study
DOI <- "10.1016/j.lwt.2021.110877"

# information on the platform and arrangement

platform <- "Illumina" # (or set to "Illumina" or "Ion_Torrent" or "F454")
paired_end <- T # set to true for paired end, false for single end
if (!paired_end) overlapping <- T # needed to run species assignment for SILVA

# Frequently used primer sets: --------------------------------------------

# V1-V3 
# Gray28F (5′-TTTGATCNTGGCTCAG) 16 bp and Gray519r (5′-GTNTTACNGCGGCKGCTG) 18 bp, 509 bp
# 27F 5′ AGAGTTTGATCMTGGCTCAG3’ — 519R 5′ GWATTACCGCGGCKGCTG3′
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
# Primers: --------------------------------------------
# Pro341F:  5′-CCTACGGGNBGCASCAG -3′  
# Pro805R:  5′-GACTACNVGGGTATCTAATCC -3′  

# expected amplicon length 469 including primers
primer_f <- "338F" # ACTCCTACGGGAGGCAGCAG
primer_r <- "806R"  # ACTCCTACGGGAGGCAGCAG
target1 <- "16S_DNA"
target2 <- region

# V3-V4
# forward primer 17 bp
# reverse primer 20 bp
# expected length 469

# create sub-directory ----------------------------------------------------------
# sub-directory data must already be in the wd

fastq_path <- file.path("data", "fastq")
filt_path <- file.path("data", "filtered")

# if(!file_test("-d", fastq_path)) {
#  dir.create(fastq_path)
# only needed if you want to create the directory and move the data from somewhere else
file_list <- list.files("./data/fastq")
file_list <- list.files(fastq_path)
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
sample.names <- dplyr::case_when(
  all(str_detect(basename(fnFs),"_")) ~ sapply(strsplit(basename(fnFs), "_", fixed = T), `[`, 1),
  all(str_detect(basename(fnFs),".")) ~ sapply(strsplit(basename(fnFs), ".", fixed = T), `[`, 1),
)

# use accession list to control which sequences will be processed
# only sequences in the accn_list file will be processed
use_accn_list <- F
accn_list_name <- "SraAccList_16S.txt" # need to adapt this to your own file
if(use_accn_list) {
  accn_list <- pull(read_tsv(accn_list_name, col_names = F),1)
  seq_to_process <- sample.names %in% accn_list
  fnFs <- fnFs[seq_to_process]
  if(paired_end) fnRs <- fnRs[seq_to_process]
  sample.names <- sample.names[seq_to_process]
}

if(keep_time) tic("\nReading sequences")
# check for occurrence of primers and adapters on a sample of forward and
# reverse sequences
# get a sample of 6 sequences
sampleFs <- if(length(fnFs)>=6) sample(fnFs,6) else fnFs[1:length(fnFs)]
# works with .fastq and fastq.gz
myFwsample <- ShortRead::readFastq(sampleFs)
head(sread(myFwsample),10)
tail(sread(myFwsample),10)
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

if(play_audio) beep(sound = 6)
if(keep_time) toc()

# If you have PacBio data you should consider using dada2::removePrimers()
# if you want to see a report browseURL(report(qa(sampleFs))) (or replace with fnFs and fnRs)
# it may be useful if there is a mismatch in reads between forwards and reverse
# but it is slow


rm(myFwsample)
if (paired_end) rm(myRvsample)

# creates quality profile plots --------------------------------

if(keep_time) tic("\ncreate quality profile plots")

# adapt this if you want to pick specific runs
plot_x_limit <- round(ave_seq_length/50)*50+25
toplot_fwd <- sampleFs
qplotfwd <- plotQualityProfile(toplot_fwd) +
  scale_x_continuous(limits = c(0,plot_x_limit), 
                     breaks = seq(0,plot_x_limit,25)) + 
  ggtitle("Fwd") + 
  theme(panel.grid.major.y = element_line(colour="grey75", size=0.5, linetype = 3),
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
    theme(panel.grid.major.y = element_line(colour="grey75", size=0.5, linetype = 3),
          axis.text.x = element_text(angle = 90, hjust = 1))
  qplotrev
}

# need to remove primers+6 bp; forward very good quality up to 295 bp; reverse up to 250

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
if(play_audio) beep(sound = 6)


# filtering and trimming --------------------------------------------------

# creates directory data/filtered if it does not exist
if(!file_test("-d", filt_path)) dir.create(filt_path)
# creates file paths for filtered sequences
filtFs <- file.path(filt_path, basename(fnFs))

if(paired_end) filtRs <- file.path(filt_path, basename(fnRs))


# truncf position at which forward sequences will be truncated
# truncr idem for reverse
# quality score 30 1 error in 100, 40 1 in 10000, 20 1 in 100
# forward need to discard the first 20 bases and cut at 450 (there will be a lot of loss)
# reverse not applicable here
truncf<- 290
truncr<- 275 # NULL if not paired end
trim_left = c(25,26) # use a length 2 vector c(x,y) if paired end or a single number if not
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
minfracloss <- min((out[,1]-out[,2])/out[,2], na.rm = TRUE)
maxfracloss <- max((out[,1]-out[,2])/out[,2], na.rm = TRUE)
medfracloss <- median((out[,1]-out[,2])/out[,2], na.rm = TRUE)

if(play_audio) beep(sound = 6)

write_tsv(filter_and_trim_par, "filtertrimpars_conc.txt")

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


# learning error rates ----------------------------------------------------

if(keep_time) tic("\nlearning error rates")

# as a default uses the first 1M sequences; if you want more nbases = xx
# another option could be randomize = T, to randomly pick samples for error
# estimates
# this may take VERY LONG
errF <- learnErrors(filtFs, multithread=TRUE, randomize = F)
if(paired_end) errR <- learnErrors(filtRs, multithread=TRUE, randomize = F)

# plot error rates
plotErrors(errF, nominalQ=TRUE) + ggtitle("Fwd")
if(paired_end) plotErrors(errR, nominalQ=TRUE) + ggtitle("Rev")
# look at how the black line follows the points

if(play_audio) beep(sound = 6)
if(keep_time) toc()

# dereplicate, infer, merge -------------------------------------------------------------

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
if(play_audio) beep(sound = 6)

# merge sequences

# If paired end, when possible merge together the inferred forward and reverse sequences.
# with a good overlap (30 bp)
# set verbose = F if you want less output

if(paired_end){
  if(keep_time) tic("\nmerging paired ends")
  minO = 25 # should be 20-30
  maxM = 0 # should be 0
  overlapping <- T # set this to true if you think there is good overlap
  # when  fwd and rev are not able to cover the full amplicon length
  # with a given overlap
  # justConcatenate = T adds 10 N between fwd and rev
  if (overlapping){
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, 
                          maxMismatch = maxM,
                          minOverlap = minO,
                          verbose=TRUE)
  } else {
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,
                          justConcatenate = T,
                          verbose=TRUE)
  }
  
  if(keep_time) toc()
  
} else {
  mergers <- dadaFs
}

if(play_audio) beep(sound = 6)

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
barplot(lengthdistr, main = "sequence length distirbution")

if(keep_time) toc()


# filter by size ----------------------------------------------------------

# to remove sequences much shorter or longer than expected (approx 420-430 for 
# V3-V4, 268-262 for V4, 388-392 for V4-V5 and 470-475 for V1-V3; shorter 
# sequences are usually mitochondria). Needs to be inspected.

seqtab_f <- seqtab.all[,nchar(colnames(seqtab.all)) %in% seq(429, 431)]

filt_seqs <- sum(seqtab_f)/sum(seqtab.all) 

# look at the abundance distribution
qplot(log10(colSums(seqtab_f)), main = "sequence length distribution", 
      xlab = "log10(seq_count)", ylab = "number of individual sequences")
# several singletons, will probably go away during chimera removal

save.image(file = str_c(Study,".Rdata"))

# remove bimera -------------------------------------------------------------
if(keep_time) tic("\nremove bimeras")
# may be very slow for large studies, you may want to check with a shorter version
# seqtab.nochim <- removeBimeraDenovo(seqtab_f[,1:10000])
cat("\nsequences prior to bimera removal\n")
dim(seqtab_f)
seqtab.nochim <- removeBimeraDenovo(seqtab_f, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
beep(sound=6)
# better save workspace here
save.image(file = str_c(Study,".Rdata"))
# recheck length distribution and see if singletons are still there
qplot(log10(colSums(seqtab.nochim)), main = "sequence distr., post bimera removal", xlab = "log10(seqn)")
if(keep_time) toc()


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
cat("The sequences after bimera removal are: ", dim(seqtab.nochim), "\n")

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

if(keep_time) tic("\nassign taxonomy, genus")

taxdb_dir <- "../tax_db" # change this if the tax databases are elsewhere
list.files(taxdb_dir)


# assignment with SILVA

# SILVA
ref_fasta <- paste(taxdb_dir, "/silva_nr99_v138_1_train_set.fa", sep="")
taxtab <- assignTaxonomy(seqtab.nochim, refFasta = ref_fasta, multithread = TRUE)
# optionally, if the sequences do not get identified, add the option tryRC=T
# which also tries the reverse complement
if(keep_time) toc()

# do species assignment, DOES NOT WORK WITH JUST CONCATENATE

if(!paired_end | overlapping){
  if(keep_time) tic("\nassign taxonomy, species")
  sp_ass_SILVA <- paste(taxdb_dir, "/silva_species_assignment_v138_1.fa", sep="")
  taxtab <- addSpecies(taxtab, sp_ass_SILVA)
  # optionally, if the sequences do not get identified, add the option tryRC=T
  # which also tries the reverse complement
  if(keep_time) toc()
}

if(play_audio) beep(sound=6)

# create a column with row names 

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


# as an alternative you can use DECIPHER::IdTaxa, which has been reported to
# be faster and more accurate than the naïve Bayesian classifier
# NEED TO CHECK IF THIS WORKS; besides they still have the 2019 version of SILVA v138

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

# taxa can be cleaned with tidyr::separate

# not run, unnecessary
# colnames(taxtab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# save a smaller version of workspace which will be needed later
save(taxtab, seqtab.nochim, track2, Study, target, region, seq_accn, DOI,
     primer_f, primer_r, target1, target2, taxtab2, 
     file = str_c(Study,"_small.Rdata"))
if (play_audio) beep(sound = 6)


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

# Build phylogenetic tree ---------------------------------------------

# decipher can aligh only ~46k seqs
# https://www.bioconductor.org/packages/3.7/bioc/vignettes/DECIPHER/inst/doc/ArtOfAlignmentInR.pdf

# if this is >46k will crash
# but it is already prohibitive with 20k
# and usually not worth the effort

dotree <- T # avoid if overlapping == F
if(dim(seqtab.nochim)[2]>10000 | overlapping == F) {dotree <- F}

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
beep(sound = 6)
save.image(file = str_c(Study,"_small.Rdata"))

# Combine data into a phyloseq object -------------------------------------

# path to the metadata file (needs to be adapted)
metadata_path <- file.path("data", "metadata", "SraRunTable.txt.csv") 
# alternatives when using data downloades from NCBI SRA are:
# "data/metadata/SraRunInfo.txt" 
# "data/metadata/SraRunTable.txt.csv"
# "data/metadata/SraRunTable.txt"
samdf <- read_tsv(metadata_path) 
if(ncol(samdf)==1) samdf <- read_csv(metadata_path)
if(use_accn_list) samdf <- dplyr::filter(samdf, Run %in% sample.names)

# check if any cols in seqtab have 0 sums
seq_sums <- rowSums(seqtab.nochim)
runs_to_keep <- which(seq_sums>0)
runs <- rownames(seqtab.nochim)[runs_to_keep]
if(length(runs_to_keep) < dim(seqtab.nochim)[1]) {
  seqtab.nochim <- seqtab.nochim[runs_to_keep,]
  samdf <- dplyr::filter(samdf, Run %in% runs)
}

# sanity check
if(all(rownames(seqtab.nochim) %in% samdf$Run)){
  cat("\nsamples in sequence table match samples in metadata\n")
} else {
  cat("\nWARNING samples in sequence table DO NOT match samples in metadata\n")
}
# rownames(seqtab.nochim) <- sapply(strsplit(rownames(seqtab.nochim), "_", fixed = T), `[`, 1)

# add final number of sequences and number of issues from track2
samdf <- left_join(samdf, select(track2, Run = sample_name, n_reads2 = nonchim, 
                                 n_issues))

samdf <- as.data.frame(samdf)
rownames(samdf) <- samdf$Run

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
# extract the average read length as an integer
read_length <- round(mean(nchar(rownames(ttab)), na.rm = T)) 

# adapt this or sett manually a comma delimited string of countries.
loc_list <- str_c(flatten_chr(distinct(samdf, geo_loc_name_country)) ,sep =",")

# put together and save study info
study <- tibble(target = target, region = region, platform = instrument,
                read_length_bp = read_length, seq_center = seq_center,
                tax_database = "SILVA v138_1", Seq_accn = seq_accn,
                samples = n_samples, DOI = DOI, geoloc = loc_list, 
                primer_f, primer_r, overlapping, paired_end)
# saves study info
write_tsv(study, str_c(Study,"_study.txt"))


# format sample information for FMBN -------------------------------------------

# this needs to be adapted for each study
# check naming of the geoloc info

samples <- samples %>%
  mutate(description = str_c(Isolation_source, Sample.Name, sep =", "))
samples <- samples %>%
  mutate(Sample_Name = Run) 


# create label2 (to avoid numbers as first char.; s. can be removed later with
# tidyr::separate)
samples <- samples %>%
  mutate(label2 = str_c("s.",Sample_Name), target1 = target1, 
         target2 = target2) %>%
  select(label2, n_reads2, n_issues, description, target1, target2, 
         biosample = BioSample, SRA_Sample = BioSample, SRA_run = Run, 
         geo_loc_country = geo_loc_name_country, 
         geo_loc_continent = geo_loc_name_country_continent, lat_lon)


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
  exists("paired_end")
  
)


mylist <- list(Study_accn = Study,
               study_df = study,
               overlap = overlapping,
               pend = paired_end,
               physeq = myphseq)




if(all(check_list)){
  cat("\nAll needed objects available and ready to process\n")
  saveRDS(mylist, file = str_c(Study, "_mindata.RDS"))
} else {
  cat("One or more of the objects you need is missing, check your data before proceeding\n")
}


cat("\nSaved data for ",Study,"\n")


# Credits and copyright ---------------------------------------------------

# Most of the script is taken from https://benjjneb.github.io/dada2/tutorial.html
# with some changes and adaptations

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
