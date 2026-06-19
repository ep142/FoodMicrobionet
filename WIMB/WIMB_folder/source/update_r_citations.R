# a function for creating and/or updating the bibliography
# for the R packages used in the script; should be run in the set-up block
update_r_citations <- function(packages, file_name = "R_citations.bib", quiet = TRUE) {
  
  log_msg <- function(msg, is_warning = FALSE) {
    if (!quiet) {
      if (is_warning) warning(msg, call. = FALSE) else message(msg)
    }
  }
  
  if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
  
  bib_entries <- list()
  
  # Strip system defaults to prevent overlapping key locks
  target_packages <- unique(packages[!(packages %in% c("stats", "tools", "utils", "methods"))])
  
  for (pkg in target_packages) {
    if (system.file(package = pkg) == "") next
    
    pkg_cit <- citation(pkg)
    
    # Squeeze down the object to ONLY its first citation entry 
    if (inherits(pkg_cit, "bibentry")) {
      pkg_cit <- pkg_cit[1]
    } else if (is.list(pkg_cit)) {
      pkg_cit <- pkg_cit[[1]]
    }
    
    p_text <- as.character(toBibtex(pkg_cit))
    
    # Forcefully match everything up to the first comma and replace with the clean key format
    p_text[1] <- stringr::str_replace(p_text[1], "^@([a-zA-Z]+)\\{[^,]*", paste0("@\\1{R-", pkg))
    
    # FIXED SAFETY CHECK: Check for a trailing comma reliably using str_detect
    if (!stringr::str_detect(p_text[1], ",\\s*$")) {
      p_text[1] <- paste0(p_text[1], ",")
    }
    
    bib_entries <- c(bib_entries, list(p_text))
  }
  
  new_bib_lines <- unlist(bib_entries)
  
  log_msg(paste0("Writing clean bibliography: ", file_name))
  writeLines(new_bib_lines, file_name)
  
  return(invisible(file_name))
}