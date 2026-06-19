# get_FMBN
# a function for loading FMBN prior to data analysis
# the function:
# 1. uses tools::R_user_dir() to find a single, global cache folder across different computers
# 2. uses use_local (default TRUE) to decide if a cached local copy should be used
# 3. uses force_download (default FALSE) to bypass the cache and download fresh
# 4. uses cache_max_days (default 180) to check staleness of cached copy
# 5. uses do_cache (default TRUE) to save downloaded copies locally for next time
# 6. guards against contradictory arguments (cannot have use_local=TRUE and force_download=TRUE)
# NOTE: reproducibility is guaranteed by pinning FMBN_path to a specific commit hash,
# not by cache settings.

# USAGE:
# FMBN_plus <- get_FMBN()  # Always assign the return value to a variable

# v1.5 2026/06/06
get_FMBN <- function(
    FMBN_repo = "ep142/FoodMicrobionet/",  
    FMBN_path = "29db4568f9dbcbc1ebc39a44b7b538fd708f4950/the_real_thing/FMBN_5_1/FMBN_plus.rds",
    use_local = TRUE,
    force_download = FALSE,
    cache_max_days = 180,
    do_cache = TRUE,
    local_FMBN_path = NULL # Determined dynamically below based on the active OS
) {
  
  # 1. If the user didn't specify a path, find the OS-specific global cache folder
  # tools::R_user_dir is part of base R and works natively on Windows, macOS, and Linux
  if (is.null(local_FMBN_path)) {
    global_cache_dir <- tools::R_user_dir("FoodMicrobionet", which = "cache")
    local_FMBN_path  <- file.path(global_cache_dir, "FMBN_plus.rds")
  }
  
  # Guard against contradictory arguments
  if (use_local && force_download) {
    stop(
      "Contradictory arguments: cannot have use_local = TRUE and force_download = TRUE.\n",
      "Set use_local = FALSE if you want to force a download.",
      call. = FALSE
    )
  }
  
  # Decide: do we attempt to load from cache?
  load_from_cache <- use_local && file.exists(local_FMBN_path) && !force_download
  
  # If attempting cache, check staleness
  if (load_from_cache) {
    file_age_days <- as.numeric(difftime(Sys.time(), file.info(local_FMBN_path)$mtime, units = "days"))
    if (file_age_days >= cache_max_days) {
      message("Cache is stale (", round(file_age_days, 1), " days >= ", 
              cache_max_days, " days). Downloading fresh copy...")
      load_from_cache <- FALSE
    } else {
      message("Loading cached FMBN from: ", local_FMBN_path, 
              " (", round(file_age_days, 1), " days old)")
    }
  } else if (use_local && !file.exists(local_FMBN_path)) {
    message("No cache found at: ", local_FMBN_path, ". Downloading...")
  }
  
  if (load_from_cache) {
    FMBN_plus <- readRDS(local_FMBN_path)
  } else {
    # Download from GitHub
    FMBN_raw_URL <- paste0("https://raw.githubusercontent.com/", FMBN_repo, FMBN_path)
    message("Downloading FMBN from GitHub...")
    
    con <- url(FMBN_raw_URL, method = "libcurl")
    on.exit(close(con), add = TRUE)
    FMBN_plus <- readRDS(con)
    
    # Optionally save to cache
    if (do_cache) {
      if (!dir.exists(dirname(local_FMBN_path))) {
        message("Creating a new global cache directory at: ", dirname(local_FMBN_path))
        dir.create(dirname(local_FMBN_path), showWarnings = FALSE, recursive = TRUE)
      }
      saveRDS(FMBN_plus, local_FMBN_path)
      message("Saved fresh FMBN copy to cache: ", local_FMBN_path)
    }
  }
  
  return(FMBN_plus)
}
