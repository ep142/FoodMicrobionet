#' Load FoodMicrobionet (FMBN) from a global cache or GitHub
#'
#' Loads the FMBN_plus object used in any of the FoodMicrobionet analysis
#' pipelines. By default the function checks a single, OS-appropriate 
#' global cache folder (found with [tools::R_user_dir()], so
#' the same cache is reused across projects and RStudio sessions on a given
#' machine) and only re-downloads from GitHub when there is no cached copy,
#' the cache is stale, or a fresh download is explicitly requested.
#'
#' @details
#' **Reproducibility note.** `FMBN_path` pins a specific GitHub commit SHA,
#' so a given default value always points at the same immutable snapshot of
#' FMBN. This is deliberate: once an analysis is finalized you should set
#' `force_download = FALSE` and increase `cache_max_days` to a large value,
#' so the exact same file is reloaded on every re-run. However, this also
#' means the pinned SHA must be **bumped manually** whenever a new FMBN
#' release should become the default; forgetting to do so will silently keep
#' serving an outdated version to anyone calling `get_FMBN()` without
#' overriding `FMBN_path`.
#'
#' **Memory note.** `FMBN_plus` (version 5.1.x) is roughly 13 MB compressed
#' on disk but expands to roughly 80-90 MB in memory once loaded (mostly
#' the `edges_B`/`edges_F` tables, which together hold well over a million
#' rows). That is a manageable footprint on its own, but it is the base load
#' for every downstream step, so it is worth keeping in mind when budgeting
#' RAM for the rest of a pipeline run on constrained machines (e.g. 8 GB
#' RAM). The function reports the loaded object's actual size on each call
#' so this can be checked empirically rather than assumed.
#'
#' @param FMBN_repo Character. GitHub repository path (owner/repo/), passed
#'   through to the raw-content URL. Default `"ep142/FoodMicrobionet/"`.
#' @param FMBN_path Character. Path within the repo to the `.rds` file,
#'   normally including a pinned commit SHA for reproducibility (see
#'   Details).
#' @param use_local Logical. If `TRUE` (default), attempt to use a cached
#'   local copy before downloading. Cannot be `TRUE` at the same time as
#'   `force_download`.
#' @param force_download Logical. If `TRUE`, bypass the cache entirely and
#'   download a fresh copy. Default `FALSE`.
#' @param cache_max_days Numeric. Maximum age, in days, of a cached copy
#'   before it is considered stale and a fresh download is triggered.
#'   Default `180`.
#' @param do_cache Logical. If `TRUE` (default), save a freshly downloaded
#'   copy to the local cache for next time.
#' @param local_FMBN_path Character or `NULL`. Path to the local cache file.
#'   If `NULL` (default), determined automatically from
#'   `tools::R_user_dir("FoodMicrobionet", which = "cache")`, giving a single
#'   global cache shared across projects.
#'
#' @return The `FMBN_plus` object (a list; see FoodMicrobionet documentation
#'   for its internal structure).
#'
#' @examples
#' \dontrun{
#' # Typical use: load from cache if available and not stale
#' FMBN_plus <- get_FMBN()
#'
#' # Force a fresh download without caching it (e.g. to check for a new release)
#' FMBN_plus <- get_FMBN(use_local = FALSE, force_download = TRUE, do_cache = FALSE)
#'
#' # Finalize an analysis: pin to the cached copy indefinitely
#' FMBN_plus <- get_FMBN(force_download = FALSE, cache_max_days = 36500)
#' }
#'
#' @note Documentation drafted with LLM assistance and reviewed/edited by
#'   the author.
#'
#' v1.5 2026/07/25
#' @export
get_FMBN <- function(
    FMBN_repo = "ep142/FoodMicrobionet/",
    FMBN_path = "29db4568f9dbcbc1ebc39a44b7b538fd708f4950/the_real_thing/FMBN_5_1/FMBN_plus.rds",
    use_local = TRUE,
    force_download = FALSE,
    cache_max_days = 180,
    do_cache = TRUE,
    local_FMBN_path = NULL # Determined dynamically below based on the active OS
) {

  # 0. Check required namespaces up front; fail fast with a clear message
  #    rather than silently installing packages on the user's system.
  missing_pkgs <- c("rlang", "stringr")[
    !vapply(c("rlang", "stringr"), requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_pkgs) > 0) {
    stop(
      "get_FMBN() requires the following package(s) which are not installed: ",
      paste(missing_pkgs, collapse = ", "), ".\n",
      "Install with: install.packages(c(", 
      paste0('"', missing_pkgs, '"', collapse = ", "), "))",
      call. = FALSE
    )
  }

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
      rlang::inform(
        paste0("Cache is stale (", round(file_age_days, 1), " days >= ",
               cache_max_days, " days). Downloading fresh copy...")
      )
      load_from_cache <- FALSE
    } else {
      rlang::inform(
        paste0("Loading cached FMBN from: ", local_FMBN_path, " (", round(file_age_days, 1), " days old)")
      )
    }
  }

  if (load_from_cache) {
    FMBN_plus <- tryCatch(
      readRDS(local_FMBN_path),
      error = function(e) {
        stop(
          "Failed to read cached FMBN copy at ", local_FMBN_path, ".\n",
          "The file may be corrupted (e.g. from an interrupted previous download).\n",
          "Try again with force_download = TRUE, use_local = FALSE to bypass it.\n",
          "Original error: ", conditionMessage(e),
          call. = FALSE
        )
      }
    )
  } else {
    # Download from GitHub
    FMBN_raw_URL <- stringr::str_c("https://raw.githubusercontent.com/", FMBN_repo, FMBN_path)
    rlang::inform("Downloading FMBN from GitHub...")

    FMBN_plus <- tryCatch({
      con <- url(FMBN_raw_URL, method = "libcurl")
      on.exit(close(con), add = TRUE)
      readRDS(con)
    }, error = function(e) {
      stop(
        "Failed to download FMBN from:\n  ", FMBN_raw_URL, "\n",
        "Check your internet connection, and if the problem persists, check whether\n",
        "FMBN_repo/FMBN_path are still valid (a pinned commit SHA can go stale if the\n",
        "upstream repository history changes).\n",
        "Original error: ", conditionMessage(e),
        call. = FALSE
      )
    })

    # Optionally save to cache
    if (do_cache) {
      # Check if the folder doesn't exist yet before creating it
      if (!dir.exists(dirname(local_FMBN_path))) {
        rlang::inform(
          paste0("Creating a new global cache directory at: ", dirname(local_FMBN_path))
        )
        dir.create(dirname(local_FMBN_path), showWarnings = FALSE, recursive = TRUE)
      }

      # Save the file to the central directory
      saveRDS(FMBN_plus, local_FMBN_path)
      rlang::inform(
        paste0("Saved fresh FMBN copy to cache: ", local_FMBN_path)
      )
    }
  }

  # Sanity check before handing the object back: confirm the core top-level
  # tables are present, so a corrupted/partial file or an unexpected format
  # fails fast here with a clear message, rather than surfacing later as a
  # cryptic error deep in some downstream function. This checks presence of
  # the tables the pipeline actually depends on; it deliberately does not
  # pin exact row/column counts, since those grow with each FMBN release.
  required_elements <- c(
    "studies", "primers", "samples_B", "samples_F",
    "edges_B", "edges_F", "taxa", "version"
  )
  missing_elements <- setdiff(required_elements, names(FMBN_plus))
  if (!is.list(FMBN_plus) || length(missing_elements) > 0) {
    stop(
      "The loaded FMBN object does not look like a valid FMBN_plus object.\n",
      if (length(missing_elements) > 0) {
        paste0("Missing expected element(s): ", paste(missing_elements, collapse = ", "), ".\n")
      },
      "The cached or downloaded file may be corrupted, from an incompatible ",
      "FMBN version, or in an unexpected format.",
      call. = FALSE
    )
  }

  rlang::inform(
    paste0(
      "Loaded FMBN_plus (", FMBN_plus$version, "), ",
      format(object.size(FMBN_plus), units = "MB", standard = "SI"), " in memory"
    )
  )

  return(FMBN_plus)
}

# ---------------------------------------------------------------------------
# Revision note (LLM-assisted, reviewed and edited by the author):
# Claude (Anthropic) contributed the following changes to v1.4, verified
# against the actual FMBN_plus.rds file where noted:
#   - Replaced silent install.packages() with an upfront check that fails
#     fast, listing any missing required packages (rlang, stringr).
#   - Wrapped cache-read and download in tryCatch() for clearer error
#     messages, and used on.exit() to ensure the download connection is
#     closed even if readRDS() fails mid-read.
#   - Added a post-load sanity check on FMBN_plus's top-level structure,
#     tested against the real object and against a deliberately malformed
#     one (see conversation of 2026/07/25).
#   - Added an in-memory size + version report on load; caught and fixed a
#     non-existent utils::object_size() call in the process, corrected to
#     base R's object.size().
#   - Added Roxygen documentation, including empirically measured on-disk
#     vs. in-memory size (~13 MB vs. ~83 MB for v5.1.2) and a note on the
#     manual-SHA-bump reproducibility trap.
# 2026/07/25
# ---------------------------------------------------------------------------
