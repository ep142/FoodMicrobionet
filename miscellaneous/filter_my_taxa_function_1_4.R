# filter_my_taxa_1_4 13/02/26
# filter_my_taxa ----------------------------------------------------------

# a function for filtering taxa based on prevalence and abundance thresholds or
# permutation filtering (see https://doi.org/10.1093/biostatistics/kxy020)
# takes a phyloseq object and other options and returns a list
# _myphyseq_ is the phyloseq object to process
# _name_ is the physeq object name (will be taken from the list during processing)
# _prevfilter_ is either "prevab" (for a prevalence and abundance filter)
# or "PERFect" for a permutation filter;
# if prevfilter is PERFect the function will try a PERFect filter first,
# if this fails, it will try the prevab filter with default options.
# if this too fails will return a message and the original unfiltered phyloseq
# _PERFect_options_ are the options for the permutation filter, can be passed as list
# _exclude_0_ will exclude 0 in calculating mean and median abundance, default is F
# _prevthreshold_ the prevalence threshold, as fraction of samples, default 0.02
# _passboth_ a flag, T if both the prevalence and abundance filter should be passed
# _abthreshold_ the abundance threshold, default 0.001
# _filenm_ a prefix which will be used in filenames, default "prevfilter"
# _colvar_pattern_ is the level of taxonomic aggregation in the prevalence
# and abundance plot, can be "[P|p]hylum" (the default) or "[C|c]lass",
# _gres_ the resolution of graphs, default is 150 dpi
# _gtype_ graph type ("jpeg", "tiff", depending on your system), default is jpeg
# _saveplot_ a flag: should the prevalence and abundance plot be saved (default = T)
# _printplot_ a flag: should the prevalence and abundance plot be printed (default = T)
# _outfolder_ the output folder, default is the working directory
# _savepat_ a flag: should the prevalence and abundance table be saved (default = T)

# the function works with defaults and returns a list

# install/load packages 
.scran_packages <- c("tidyverse", "tictoc", "crayon")
.sbioc_packages <- c("BiocManager", "phyloseq")
.attached_pckgs <- .packages()

# install and load Bioconductor packages if necessary
if(!all(.sbioc_packages %in% .attached_pckgs)){
  .sinst <- .sbioc_packages %in% installed.packages()
  if(any(!.sinst)) {
    if(!.sinst[1]) {
      install.packages("BiocManager")
      .sinst <- .sbioc_packages %in% installed.packages()
    }
    if(any(!.sinst[2:length(.sinst)])) {
      BiocManager::install(.sbioc_packages[!.sinst], ask = F)
    }
  }
  sapply(.sbioc_packages, require, character.only = TRUE)
}

# install and load CRAN packages if necessary 
if(!all(.scran_packages %in% .attached_pckgs)){
  .sinst <- .scran_packages %in% installed.packages()
  if(any(!.sinst)) {
    install.packages(.scran_packages[!.sinst])
  }
  sapply(.scran_packages, require, character.only = TRUE)
}

# PERFect is not available any more from Bioconductor
# a file with PERFect functions must be available in the source folder
if(file.exists(file.path("source", "PERFect_functions.R"))) {
  source(file.path("source", "PERFect_functions.R"))
} else {
  stop(paste0(
    "The source file 'PERFect_functions.R' was not found in the Source folder."
  ))
}

filter_my_taxa <- function(myphyseq, 
                           verbose_output = T,
                           name = "phyloseqobj",
                           prevfilter = "prevab",
                           PERFect_options = list(
                             method = "sim", # "sim" or "perm", perm takes longer
                             perf_alpha = 0.1, # this is the default
                             perf_k = 10000, # this is the number of permutation, default value
                             rollm = F, # this is the rolling mean argument
                             perf_alg = "fast", # algorithm for permutation filtering fast or full, full is slower
                             perf_hist = FALSE, # the hist option in PERFect_perm
                             vrb_out = TRUE # if extra output needs to be produced
                           ),
                           exclude_0 = F,
                           prevthreshold = 0.02,
                           passboth = T,
                           abthreshold = 0.001,
                           filenm = "prevfilter",
                           colvar_pattern = "[P|p]hylum",
                           gres = 150,
                           gtype = "jpeg",
                           saveplot = T,
                           printplot = T,
                           outfolder = ".",
                           savepat = T){
  
  # check the object provided for filtering is a phyloseq
  if(class(myphyseq)!="phyloseq") {
    errorm <- paste0(
      "\nThe class of the object you provided is ",
      class(myphyseq),
      "\nYou must provide a phyloseq object!\n")
    stop(
      errorm,
      call. = FALSE
    )
  }
  
  # check if the output directory exists, if not creates it, to avoid errors
  if(!outfolder == "."){
    if(!dir.exists(outfolder)) dir.create(outfolder)
  }
  
  # creates list for the results, the object which will be returned by the function
  tax_filter_results <- list(phyobj = NULL,
                             phyobjname = NULL,
                             prevabplot = NULL,
                             prevabtable = NULL,
                             filter_type = NULL,
                             exclude_0 = exclude_0)

  # check if any taxon has 0 counts and removes it
  if(any(taxa_sums(myphyseq)==0)){
    ntaxa_pre <- ntaxa(myphyseq)
    myphyseq <- filter_taxa(myphyseq, function(x) sum(x)>0, prune = T)
    ntaxapost <- ntaxa(myphyseq)
    taxaremoved <- ntaxa_pre-ntaxapost
    if(verbose_output) red(cat("\n",taxaremoved," taxa were removed because of 0 sums","\n"))
  }
  
  # check if any samples has 0 counts and removes it
  if(any(sample_sums(myphyseq)==0)){
    nsamples_pre <- nsamples(myphyseq)
    myphyseq <- subset_samples(myphyseq, sample_sums(myphyseq)>0)
    nsamplespost <- nsamples(myphyseq)
    samplesremoved <- nsamplespre-nsamplespost
    if(verbose_output) red(cat("\n",samplesremoved," samples were removed because of 0 sums","\n"))
  }
  
  # get the OTU table
  # get the option on the OTU table
  taxa_are_rows_flag <- myphyseq@otu_table@taxa_are_rows
  # obtain a prevalence and abundance plot
  OTUmatrixf <- as(otu_table(myphyseq), "matrix")
  if (taxa_are_rows_flag) {
    OTUmatrixf <- t(OTUmatrixf)
  }
  
  OTUmatrixf_relab <- OTUmatrixf / rowSums(OTUmatrixf)
  # calculate the prevalence on all columns
  prevdf <- apply(
    X = OTUmatrixf,
    MARGIN = 2,
    FUN = function(x) {
      sum(x > 0)
    }
  )
  # calculate minimum relative abundance
  min_rel_ab <- apply(
    X = OTUmatrixf_relab,
    MARGIN = 2,
    FUN = function(x) {
      min(x)
    }
  )
  # calculate max. rel abundance
  max_rel_ab <- apply(
    X = OTUmatrixf_relab,
    MARGIN = 2,
    FUN = function(x) {
      max(x)
    }
  )
  # replace 0 with NA if exclude_0 is T
  if(exclude_0){
    OTUmatrixf_relab[OTUmatrixf_relab == 0] <- NA_real_
  }
  
  # calculate median relative abundance
  med_rel_ab <- apply(
    X = OTUmatrixf_relab,
    MARGIN = 2,
    FUN = function(x) {
      median(x, na.rm = T)
    }
  )
  
  # calculate mean relative abundance
  mean_rel_ab <- apply(
    X = OTUmatrixf_relab,
    MARGIN = 2,
    FUN = function(x) {
      mean(x, na.rm = T)
    }
  )
  
  prevdf <- data.frame(
    Prevalence = prevdf,
    TotalAbundance = colSums(OTUmatrixf),
    min_rel_ab = min_rel_ab,
    max_rel_ab = max_rel_ab,
    med_rel_ab = med_rel_ab,
    mean_rel_ab = mean_rel_ab
  )
  
  # merge taxonomic information
  taxa_metadata <-
    as.data.frame(as(tax_table(myphyseq),"matrix")) |> 
    rownames_to_column(var = "label")
  tranks <- rank_names(myphyseq)
  prevdf <- prevdf |>
    rownames_to_column(var = "label") %>%
    left_join(., select(taxa_metadata, label, any_of(tranks)))
  
  # apply the filter 
  # "prevab" is a prevalence and abundance filter
  # "PERFect" is a permutation filter, see https://doi.org/10.1093/biostatistics/kxy020.
  if (prevfilter == "PERFect") {
    # use PERFect
    # create the count data frame
    if (taxa_are_rows_flag) {
      counts <- as.data.frame(t(otu_table(myphyseq)))
    } else {
      counts <- as.data.frame(otu_table(myphyseq))
    }
    if (PERFect_options$method == "sim") {
      tic("start PERFect_sim")
      res_sim <- try(PERFect_sim(X = counts, alpha = PERFect_options$perf_alpha))
      if (class(res_sim) != "try-error") {
        if (PERFect_options$vrb_out) {
          p <- pvals_Plots(
            PERFect = res_sim,
            X = counts,
            quantiles = c(0.25, 0.5, 0.8, 0.9),
            alpha = PERFect_options$perf_alpha
          )
          print(p$plot + ggtitle(str_c(
            filenm, " Simultanenous Filtering"
          )))
        }
        OTUtokeep <- colnames(res_sim$filtX)
      } else {
        OTUtokeep <- NULL
      }
      toc()
    } else {
      tic("start PERFect_perm")
      res_sim <- try(PERFect_perm(
        X = counts,
        alpha = PERFect_options$perf_alpha,
        algorithm = PERFect_options$perf_alg,
        hist = PERFect_options$perf_hist,
        rollmean = PERFect_options$rollm
      ))
      if (class(res_sim) != "try-error") {
        if (PERFect_options$vrb_out) {
          p <- pvals_Plots(
            PERFect = res_sim,
            X = counts,
            quantiles = c(0.25, 0.5, 0.8, 0.9),
            alpha = PERFect_options$perf_alpha
          )
          print(p$plot + ggtitle(
            str_c(
              filenm,
              "Permutation Filtering",
              PERFect_options$perf_alg,
              "algorithm",
              sep = " "
            )
          ))
        }
        OTUtokeep <- colnames(res_sim$filtX)
      } else {
        OTUtokeep <- NULL
      }
      
      toc()
    }
    if (is.null(OTUtokeep)) {
      errorm <- red("\nPERFect filter failed, switching to prevab\n")
      cat(errorm)
      prevfilter <- "prevab"
    } else {
      tax_filter_results$filter_type <- paste0(prevfilter, PERFect_options$method)
    }
  } 
  
  # this is executed if the prevfilter is prevab or if the PERFect filter fails
  if(prevfilter == "prevab"){
    pass_prev_filter <-
      dplyr::filter(select(prevdf, label, Prevalence),
                    Prevalence > floor(nsamples(myphyseq) *
                                         prevthreshold)) %>%
      pull(label)
    if (passboth) {
      OTUtokeep <- intersect(names(which(max_rel_ab >= abthreshold)),
                             pass_prev_filter)
    } else {
      OTUtokeep <- union(names(which(max_rel_ab >= abthreshold)),
                         pass_prev_filter)
    }
    if (is.null(OTUtokeep)) {
      errorm <- red("\nprevab filter failed, returning all taxa\n")
      cat(errorm)
      tax_filter_results$filter_type <- "No filter"
      OTUtokeep <- pull(prevdf, label)
    } else {
      tax_filter_results$filter_type <- paste0(prevfilter)
    }
  }

  # apply the filter to the prevdf table
  prevdf <- prevdf  |>
    mutate(
      relAbundance = TotalAbundance / sum(TotalAbundance),
      pass_filters = ifelse(label %in% OTUtokeep, "T", "F")
    )
  # apply the filter and put the result in the return list
  tax_filter_results$phyobj <- prune_taxa(OTUtokeep, myphyseq)
  tax_filter_results$phyobjname <- name
  
  # prevalence vs abundance plot
  OTUmatrixf <-
    OTUmatrixf[, which(colnames(OTUmatrixf) %in% OTUtokeep)]
  # fraction of remaining sequences
  
  f_seq_ret <-
    round(sum(sample_sums(tax_filter_results$phyobj)) / sum(sample_sums(myphyseq)), 4)
  
  # make a plot

  # original number of taxa
  ntaxa_prefilt <- ntaxa(myphyseq)
  subtitle_text <-
    paste(
      "using the filters you retain ",
      length(OTUtokeep),
      " taxa (triangles) out of ",
      ntaxa_prefilt,
      " (",
      f_seq_ret * 100,
      "% of init. seqs.)",
      sep = ""
    )
  # a prevalence and abundance plot change rank for color and facet as appropriate
  colorvar <- rank_names(myphyseq)[which(str_detect(rank_names(myphyseq), colvar_pattern))]
  title_text <- paste("Prevalence vs. abundance, by ", colorvar, ", ", name, sep = "")
  prev_ab_plot <-
    ggplot(prevdf,
           aes(
             x = TotalAbundance,
             y = Prevalence / nrow(OTUmatrixf),
             shape = as.factor(pass_filters),
             color = .data[[colorvar]]
           )) +
    geom_point(size = 2, alpha = 0.7) +
    facet_wrap(~ .data[[colorvar]]) +
    geom_hline(
      yintercept = ifelse(prevfilter == "prevab", prevthreshold, 0),
      alpha = 0.5,
      linetype = 2
    ) +
    labs(
      x = "total abundance",
      y = "Prevalence [Frac. Samples]",
      shape = 'pass ab. threshold',
      title = title_text,
      subtitle = subtitle_text
    ) +
    scale_x_log10() +
    scale_y_continuous(minor_breaks = seq(0, 1, 0.05)) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 90)
    )
  if(printplot) print(prev_ab_plot)
  # add the plot to the list
  tax_filter_results$prevabplot <- prev_ab_plot
  
  if(saveplot) ggsave(prev_ab_plot,
                      file = paste(file.path(outfolder,filenm), 
                                   "_prevabg.", gtype, sep=""),
                      device = gtype,
                      dpi = gres, 
                      width = 7,
                      height = 7)
  # finish up with the table
  prevdf <- prevdf |>
    mutate(relprev = Prevalence / nrow(OTUmatrixf)) %>%
    arrange(-relprev,-relAbundance)
  tax_filter_results$prevabtable <- prevdf
  
  if(savepat) {
    write_tsv(prevdf, 
              file = paste(file.path(outfolder,filenm), "_prevabt.txt",sep=""))
  }
  return(tax_filter_results)
}

# Credits and copyright notice

# The code for loading packages is derived from [this post](http://tinyurl.com/jjwyzph).  

# The code for PERFect functions is taken verbatim from 
# Smirnova E, Cao Q (2022). PERFect: Permutation filtration for microbiome data. R package version 1.12.0, https://github.com/cxquy91/PERFect.
# this is because most recent versions of Bioconductor do not support PERFect

# The code for the creation of the prevalence and abundance data frame is taken from here: 
# https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html  
# and significantly modified

# Assume that this is overall under MIT licence

# Copyright 2025, 2026 Eugenio Parente
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