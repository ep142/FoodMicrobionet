# Bioconductor pipeline for ITS: functions --------------------------------

# v1.2.3 5/3/2024

# three functions needed in the pipeline 
# the 16S pipeline only needs the first 2


# this function is only defined if you use a primer table
if(use_primer_table){
  double_check_primers <- function(primer_file = primer_table_path, 
                                   primerf_name, primerf_seq, primerr_name, primerr_seq){
    # is the primer file there?
    if(!file.exists(primer_table_path)) stop("Your primer table is not where you told me...")
    primer_table <- read_tsv(primer_table_path)
    # check the forward name
    primerf_n <- unique(pull(primer_table, primer_f_name))
    primerr_n <- unique(pull(primer_table, primer_r_name))
    if(!primerf_name %in% primerf_n) {
      warning("the name of your forward primer is not in the table you provided")
    } else {
        cat("\nThe name of your forward primer is in the table you provided\n")
        # check the sequence
        fprimer_from_table <- primer_table |>
          dplyr::filter(primer_f_name == primerf_name) |>
          pull(primer_f_seq) |>
          unique()
        if(length(fprimer_from_table)>1) stop("\nYou have more than one sequence associated to the primer name in the primer table.\n")
        if(fprimer_from_table == primerf_seq){
          cat("\nThe sequence of your forward primer matches the forward primer name\n")
        } else {
          warning("the sequence of your forward primer does not match with the table")
        }
      }
    if(!primerr_name %in% primerr_n) {
      warning("the name of your forward primer is not in the table you provided")
    } else {
      cat("\nThe name of your reverse primer is in the table you provided\n")
      # check the sequence
      rprimer_from_table <- primer_table |>
        dplyr::filter(primer_r_name == primerr_name) |>
        pull(primer_r_seq) |>
        unique()
      if(length(rprimer_from_table)>1) stop("\nYou have more than one sequence associated to the primer name in the primer table.\n")
      if(rprimer_from_table == primerr_seq){
        cat("\nThe sequence of your reverse primer matches the reverse primer name\n")
      } else {
        warning("the sequence of your reverse primer does not match with the table")
      }
    }
    if(primerf_name %in% primerf_n & primerr_name %in% primerr_n){
      primer_line_df <- dplyr::filter(primer_table, primer_f_seq == primer_f_seq & primer_r_seq == primerr_seq)
      if(nrow(primer_line_df) == 0) warning(" I cannot locate your primer pair, check your primers or  primer table (names, sequences)...")
    }
  }
}

# create different orientations of the primers
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

# mixed_merge function ----------------------------------------------------

# this function is modeled over:
# https://github.com/benjjneb/dada2/issues/537
# many of the arguments are devined or created in the master script
# save_output if true, saves some of the objects produced within the
# function for further analysis

mixed_merge <- function(
    ddFs = dadaFs, drFs = derepFs, ddRs = dadaRs, drRs = derepRs, 
    mnO = minO, save_output = T){
  # maxMismatch=0 has higher penalty for mismatch and gap for nwalgin (-64)
  cat("\n\nperfect match merging....\n")
  mergedData <- mergePairs(ddFs, drFs, ddRs, drRs, 
                           maxMismatch=0, minOverlap = mnO,
                           trimOverhang=TRUE, returnRejects=TRUE, 
                           verbose=TRUE)
  # maxMismatch > 0 has lower penalty for mismatch and gap for nwalgin (-8)
  cat("\n\nmismatch allowed merging....\n")
  mismatch_mergers <- mergePairs(ddFs, drFs, ddRs, drRs, 
                                 maxMismatch=4, minOverlap = mnO,
                                 trimOverhang=TRUE, returnRejects=TRUE, 
                                 verbose=TRUE)
  cat("\n\nconcatenation merging....\n")
  concats <- mergePairs(ddFs, drFs, ddRs, drRs, 
                        justConcatenate=TRUE, verbose=TRUE)
  
  # optionally save sample merge dfs for diff comparisions; 
  # a directory will be created if needed
  if(save_output){
    if(!dir.exists(file.path("data", "merger_dir"))){
      dir.create(file.path("data", "merger_dir"))
    }
    merger_dir <- file.path("data","merger_dir")
    sapply(names(mergedData), 
           \(x) write_tsv(mergedData[[x]], 
                          file=file.path(merger_dir, 
                                         str_c(x, "strict_mergers.txt", sep="_")))
    )
    sapply(names(mismatch_mergers), 
           \(x) write_tsv(mismatch_mergers[[x]], 
                          file=file.path(merger_dir, 
                                         str_c(x, "loose_mergers.txt", sep="_")))
    )
    # converting concat df to matrix to overcome error of - unimplemented type 'list' in 'EncodeElement' 
    sapply(names(concats), 
           \(x) write_tsv(as.data.frame(as.matrix(concats[[x]])), 
                          file=file.path(merger_dir, 
                                         str_c(x, "all_concats.txt", sep="_")))
    )
  }
  # replace mismatched or concatenated ASVs in the main mergedData
  for(i in names(mergedData)) {
    # store row index to drop certain ASVs later
    rowsToDelete = vector()
    mergedDf = mergedData[[i]]
    cat(i, "Out of total", sum(mergedDf$abundance), "paired-reads (in", nrow(mergedDf), "unique pairings), retained ")
    mismatchDf = mismatch_mergers[[i]]
    rownames(mismatchDf) = paste(mismatchDf$forward, mismatchDf$reverse, sep="_")
    concatDf = concats[[i]]
    rownames(concatDf) = paste(concatDf$forward, concatDf$reverse, sep="_")
    
    for (row in 1:nrow(mergedDf)) {
      # skipping rows that are good to go from default analysis
      if (mergedDf[row,]$accept) { next }
      
      uniquePairID = paste(mergedDf[row,]$forward, mergedDf[row,]$reverse, sep="_")
      
      # if match length is less than minO then good to go with concatenation
      if (mismatchDf[uniquePairID,]$nmatch <= mnO) {
        mergedDf[row,] = concatDf[uniquePairID,]
      }
      # ignoring ASVs with many mismatches with long overlap
      else{
        misMatchIndels = mismatchDf[row,]$nmismatch + mismatchDf[row,]$nindel
        cutOff = 0
        if ( mismatchDf[uniquePairID,]$nmatch > mnO && mismatchDf[uniquePairID,]$nmatch <= 50 ) {
          cutOff = 1
        }
        else if ( mismatchDf[uniquePairID,]$nmatch > 50 && mismatchDf[uniquePairID,]$nmatch <= 100 ) {
          cutOff = 2
        }
        else if (mismatchDf[uniquePairID,]$nmatch > 100) {
          cutOff = 3
        }
        # check if mismatches are below cut off
        if (misMatchIndels <= cutOff) {
          mergedDf[row,] = mismatchDf[uniquePairID,]
        } 
        # discard if mismatches are high in longer than 100bp match 
        # else if (cutOff==3 && misMatchIndels > 3) {
        # 	rowsToDelete = c(rowsToDelete, row)
        # }
        # concatenate if mismatches are high in overlap region remove reverse read part of the overlap region 
        else {
          trimLength = mismatchDf[uniquePairID,]$nmatch
          concatDf[uniquePairID,]$sequence = gsub(paste0("N{10}[A-z]{", trimLength, "}"), "NNNNNNNNNN", concatDf[uniquePairID,]$sequence, fixed = T)
          mergedDf[row,] = concatDf[uniquePairID,]
        }
      }
    }
    if (length(rowsToDelete) == 0){
      mergedData[[i]] = mergedDf
    } else {
      mergedData[[i]] = mergedDf[-rowsToDelete,]
    }
    cat(sum(mergedData[[i]]$abundance), "paired-reads (in", nrow(mergedData[[i]]), "unique pairings)\n")
  }
  return(mergedData)
}

# the same function but for the big data pipeline
mixed_merge_bigdata <- function(
    ddFs = ddF, drFs = derepF, ddRs = ddR, drRs = derepR, 
    mnO = minO, save_output = T){
  # maxMismatch=0 has higher penalty for mismatch and gap for nwalgin (-64)
  cat("\n\nperfect match merging....\n")
  mergedData <- mergePairs(ddFs, drFs, ddRs, drRs, 
                           maxMismatch=0, minOverlap = mnO,
                           trimOverhang=TRUE, returnRejects=TRUE, 
                           verbose=TRUE)
  # maxMismatch > 0 has lower penalty for mismatch and gap for nwalgin (-8)
  cat("\n\nmismatch allowed merging....\n")
  mismatch_mergers <- mergePairs(ddFs, drFs, ddRs, drRs, 
                                 maxMismatch=4, minOverlap = mnO,
                                 trimOverhang=TRUE, returnRejects=TRUE, 
                                 verbose=TRUE)
  cat("\n\nconcatenation merging....\n")
  concats <- mergePairs(ddFs, drFs, ddRs, drRs, 
                        justConcatenate=TRUE, verbose=TRUE)
  
  # optionally save sample merge dfs for diff comparisions; 
  # a directory will be created if needed
  if(save_output){
    if(!dir.exists(file.path("data", "merger_dir"))){
      dir.create(file.path("data", "merger_dir"))
    }
    merger_dir <- file.path("data","merger_dir")
    write_tsv(mergedData,
              file=file.path(merger_dir, str_c(sam, "strict_mergers.txt", sep="_")))
    write_tsv(mismatch_mergers,
              file=file.path(merger_dir, str_c(sam, "loose_mergers.txt", sep="_")))
    
    # converting concat df to matrix to overcome error of - unimplemented type 'list' in 'EncodeElement' 
    write_tsv(as.data.frame(as.matrix(concats)), 
              file=file.path(merger_dir, str_c(sam, "all_concats.txt", sep="_")))
  }
  # replace mismatched or concatenated ASVs in the main mergedData
  # store row index to drop certain ASVs later
  rowsToDelete = vector()
  mergedDf = mergedData
  cat("Out of total", sum(mergedDf$abundance), "paired-reads (in", nrow(mergedDf), "unique pairings), retained ")
  mismatchDf = mismatch_mergers
  rownames(mismatchDf) = paste(mismatchDf$forward, mismatchDf$reverse, sep="_")
  concatDf = concats
  rownames(concatDf) = paste(concatDf$forward, concatDf$reverse, sep="_")
  
  for (row in 1:nrow(mergedDf)) {
    # skipping rows that are good to go from default analysis
    if (mergedDf[row,]$accept) { next }
    
    uniquePairID = paste(mergedDf[row,]$forward, mergedDf[row,]$reverse, sep="_")
    
    # if match length is less than minO then good to go with concatenation
    if (mismatchDf[uniquePairID,]$nmatch <= mnO) {
      mergedDf[row,] = concatDf[uniquePairID,]
    }
    # ignoring ASVs with many mismatches with long overlap
    else{
      misMatchIndels = mismatchDf[row,]$nmismatch + mismatchDf[row,]$nindel
      cutOff = 0
      if ( mismatchDf[uniquePairID,]$nmatch > mnO && mismatchDf[uniquePairID,]$nmatch <= 50 ) {
        cutOff = 1
      }
      else if ( mismatchDf[uniquePairID,]$nmatch > 50 && mismatchDf[uniquePairID,]$nmatch <= 100 ) {
        cutOff = 2
      }
      else if (mismatchDf[uniquePairID,]$nmatch > 100) {
        cutOff = 3
      }
      # check if mismatches are below cut off
      if (misMatchIndels <= cutOff) {
        mergedDf[row,] = mismatchDf[uniquePairID,]
      } 
      # discard if mismatches are high in longer than 100bp match 
      # else if (cutOff==3 && misMatchIndels > 3) {
      # 	rowsToDelete = c(rowsToDelete, row)
      # }
      # concatenate if mismatches are high in overlap region remove reverse read part of the overlap region 
      else {
        trimLength = mismatchDf[uniquePairID,]$nmatch
        concatDf[uniquePairID,]$sequence = gsub(paste0("N{10}[A-z]{", trimLength, "}"), "NNNNNNNNNN", concatDf[uniquePairID,]$sequence, fixed = T)
        mergedDf[row,] = concatDf[uniquePairID,]
      }
    }
  }
  if (length(rowsToDelete) == 0){
    mergedData = mergedDf
  } else {
    mergedData = mergedDf[-rowsToDelete,]
  }
  cat(sum(mergedData$abundance), "paired-reads (in", nrow(mergedData), "unique pairings)\n")
  return(mergedData)
  }

