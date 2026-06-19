# A function for automatically detecting function dependencies in a
# list of functions provided by the user

get_function_dependencies <- function(functions_vector, envir = globalenv()) {
  pattern_library <- "(?:library|require|requireNamespace)\\s*\\(\\s*[\"']([a-zA-Z0-9.]+)[\"']|(?:library|require|requireNamespace)\\s*\\(\\s*([a-zA-Z0-9.]+)\\s*[,)]"
  pattern_colons  <- "([a-zA-Z0-9.]+)::"
  
  extract_all <- function(text, pattern) {
    m <- gregexpr(pattern, text, perl = TRUE)
    hits <- regmatches(text, m)
    unlist(hits)
  }
  
  detected_packages <- character(0)
  
  for (f_name in functions_vector) {
    if (!exists(f_name, envir = envir, inherits = FALSE, mode = "function")) next
    
    f_obj <- get(f_name, envir = envir)
    if (is.primitive(f_obj) || is.null(body(f_obj))) next
    
    f_code <- paste(deparse(body(f_obj)), collapse = "\n")
    
    # Extract package names from library/require calls
    lib_hits <- extract_all(f_code, pattern_library)
    lib_pkgs  <- sub(pattern_library, "\\1\\2", lib_hits, perl = TRUE)
    
    # Extract package names from pkg:: usage
    cc_hits  <- extract_all(f_code, pattern_colons)
    cc_pkgs  <- sub("::", "", cc_hits)
    
    detected_packages <- c(detected_packages, lib_pkgs, cc_pkgs)
  }
  
  sort(unique(na.omit(c("base", detected_packages[nchar(detected_packages) > 0]))))
}
