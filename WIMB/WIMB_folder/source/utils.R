# utils.R
# Shared utility functions for the Where Is My Bug? Shiny app.
# Sourced automatically by Shiny before app.R runs.

# ---------------------------------------------------------------------------
# resolve_abbr
# Given a character vector of full labels and a starting abbreviation width,
# returns a list(width, abbr) where every abbreviation is unique.
# Widens one character at a time until there are no collisions.
# ---------------------------------------------------------------------------
resolve_abbr <- function(vals, width) {
  max_len <- max(nchar(vals))
  while (width < max_len) {
    if (dplyr::n_distinct(stringr::str_trunc(vals, width, "center")) ==
        dplyr::n_distinct(vals)) break
    width <- width + 1
  }
  list(width = width,
       abbr  = stringr::str_trunc(vals, width, side = "center"))
}

# ---------------------------------------------------------------------------
# auto_detect_domain
# Given a vector of taxon names and a taxon level, looks them up in
# FMBN_plus$taxa and returns the domain ("Bacteria" or "Fungi").
# Stops with an informative message if taxa span both domains or are not found.
# ---------------------------------------------------------------------------
auto_detect_domain <- function(taxon_to_search, taxon_level, taxa_all) {
  candidate_rows <- taxa_all |>
    dplyr::filter(.data[[taxon_level]] %in% taxon_to_search)

  if (nrow(candidate_rows) == 0) {
    stop("None of the selected taxa were found in FMBN at level '",
         taxon_level, "'. Check spelling.")
  }

  detected_domains <- unique(candidate_rows$domain)

  if (length(detected_domains) > 1) {
    stop("Selected taxa span multiple domains: ",
         paste(detected_domains, collapse = " and "),
         ". All taxa must be either Bacteria or Fungi.")
  }

  detected_domains
}

# ---------------------------------------------------------------------------
# build_pooled_edges
# Core filtering and pooling pipeline, extracted from the Rmd get_data chunk.
# Returns a list with named elements: pooled_edges_sample_study, filt_samples,
# filt_studies, domain, taxon_label, nstudies_wtarget, nsamples_wtarget,
# nstudies_total, nsamples_total.
# ---------------------------------------------------------------------------
build_pooled_edges <- function(FMBN_plus, p) {

  taxa_all <- FMBN_plus$taxa

  # --- domain auto-detection ------------------------------------------------
  domain <- auto_detect_domain(p$taxon_to_search, p$taxon_level, taxa_all)

  # --- taxonomic filtering --------------------------------------------------
  my_taxa <- taxa_all |>
    dplyr::filter(domain == !!domain,
                  !is.na(.data[[p$taxon_level]]),
                  .data[[p$taxon_level]] %in% p$taxon_to_search)

  # handle NA at go_down_to level
  tax_levels    <- c("domain", "phylum", "class", "order",
                     "family", "genus", "species")
  target_level  <- p$go_down_to
  target_idx    <- match(target_level, tax_levels)
  parent_level  <- tax_levels[target_idx - 1]

  tax_level_below <- dplyr::case_when(
    p$taxon_level == "order"  ~ "family",
    p$taxon_level == "family" ~ "genus",
    p$taxon_level == "genus"  ~ "species"
  )
  if (!p$keep_NA_lowest) {
    my_taxa <- dplyr::filter(my_taxa, !is.na(.data[[tax_level_below]]))
  }

  # --- edges ----------------------------------------------------------------
  edges <- if (domain == "Bacteria") FMBN_plus$edges_B else FMBN_plus$edges_F

  if (p$remove_no_phylum) {
    no_phylum_id <- taxa_all |>
      dplyr::filter(is.na(phylum)) |>
      dplyr::pull(taxonId)
    edges <- dplyr::filter(edges, !(taxonId %in% no_phylum_id))
  }

  if (any(p$remove_chloroplast, p$remove_mitochondria, p$remove_no_phylum)) {
    edges <- edges |>
      dplyr::mutate(weight_old = weight) |>
      dplyr::group_by(sampleId) |>
      dplyr::mutate(weight = weight / sum(weight)) |>
      dplyr::ungroup()
  }

  edges_filt <- dplyr::semi_join(edges, my_taxa, by = "taxonId") |>
    dplyr::left_join(dplyr::select(my_taxa, taxonId:species), by = "taxonId")

  # fill NA at go_down_to
  if (target_idx > 1) {
    edges_filt <- edges_filt |>
      dplyr::mutate(
        !!rlang::sym(target_level) := dplyr::case_when(
          !is.na(!!rlang::sym(target_level)) ~ !!rlang::sym(target_level),
          target_level == "species"           ~ stringr::str_c(
            !!rlang::sym(parent_level), " spp."),
          TRUE                                ~ !!rlang::sym(parent_level)
        )
      )
  }

  pooled_edges <- edges_filt |>
    dplyr::group_by(sampleId, .data[[p$go_down_to]]) |>
    dplyr::summarise(weight = sum(weight), .groups = "drop") |>
    dplyr::arrange(sampleId, .data[[p$go_down_to]])

  # --- samples --------------------------------------------------------------
  filt_samples <- if (domain == "Bacteria") {
    FMBN_plus$samples_B
  } else {
    FMBN_plus$samples_F
  }

  filt_samples <- filt_samples |>
    dplyr::filter(
      !stringr::str_detect(L1, "blank"),
      !stringr::str_detect(L1, "Mock"),
      !stringr::str_detect(L1, "^Water")
    )

  if (p$s_type %in% c("Environment", "Sample")) {
    filt_samples <- dplyr::filter(filt_samples, s_type == p$s_type)
  }

  pooled_edges_sample <- dplyr::inner_join(
    pooled_edges,
    dplyr::select(filt_samples, studyId, sampleId, s_type, n_reads2,
                  n_issues, foodId, L1:target2, geo_loc_country),
    by = "sampleId"
  )

  # --- studies --------------------------------------------------------------
  filt_studies <- FMBN_plus$studies

  pooled_edges_sample_study <- dplyr::inner_join(
    pooled_edges_sample,
    dplyr::select(filt_studies, studyId, read_length_bp, Seq_accn,
                  food_group, overlapping, paired_end),
    by = "studyId"
  )

  # --- depth and NA filter --------------------------------------------------
  pooled_edges_sample_study <- pooled_edges_sample_study |>
    dplyr::filter(!is.na(.data[[p$food_grouping_variable]]),
                  n_reads2 >= p$min_n_reads)

  if (p$run_further_filters) {
    pooled_edges_sample_study <- pooled_edges_sample_study |>
      dplyr::filter(target2 %in% p$gene_targets)
  }

  # read length filter (apply only if p$domain == "Bacteria")
  if(p$domain == "Bacteria"){
    pooled_edges_sample_study <- pooled_edges_sample_study |>
    dplyr::filter(read_length_bp >= p$min_read_length)
    }

  # --- summary counts -------------------------------------------------------
  nstudies_wtarget <- dplyr::n_distinct(pooled_edges_sample_study$studyId)
  nsamples_wtarget <- dplyr::n_distinct(pooled_edges_sample_study$sampleId)
  nstudies_total   <- nrow(FMBN_plus$studies)
  nsamples_total   <- if (domain == "Bacteria") {
    nrow(FMBN_plus$samples_B)
  } else {
    nrow(FMBN_plus$samples_F)
  }

  list(
    pooled_edges_sample_study = pooled_edges_sample_study,
    filt_samples              = filt_samples,
    filt_studies              = filt_studies,
    domain                    = domain,
    taxon_label               = paste(p$taxon_to_search, collapse = ", "),
    nstudies_wtarget          = nstudies_wtarget,
    nsamples_wtarget          = nsamples_wtarget,
    nstudies_total            = nstudies_total,
    nsamples_total            = nsamples_total
  )
}

# ---------------------------------------------------------------------------
# build_prev_ab
# Computes n_samples, prev_ab_df_sel from a pooled_edges_sample_study df.
# Returns a list with n_samples and prev_ab_df_sel.
# ---------------------------------------------------------------------------
build_prev_ab <- function(pooled_edges_sample_study, filt_samples, p) {

  n_samples <- filt_samples |>
    dplyr::filter(sampleId %in% unique(pooled_edges_sample_study$sampleId)) |>
    dplyr::summarise(n = dplyr::n_distinct(sampleId),
                     .by = dplyr::all_of(p$food_grouping_variable))

  n_samples_all <- filt_samples |>
    dplyr::filter(s_type == p$s_type) |>
    dplyr::summarise(n = dplyr::n_distinct(sampleId),
                     .by = dplyr::all_of(p$food_grouping_variable)) |>
    dplyr::rename(n_all = n)

  n_samples <- dplyr::left_join(n_samples, n_samples_all,
                                 by = p$food_grouping_variable) |>
    dplyr::mutate(prev = n / n_all) |>
    dplyr::arrange(dplyr::desc(prev))

  prev_df_sel <- pooled_edges_sample_study |>
    dplyr::group_by(
      dplyr::across(dplyr::all_of(c(p$go_down_to, p$food_grouping_variable)))
    ) |>
    dplyr::summarise(abs_prev = dplyr::n_distinct(sampleId),
                     .groups = "drop") |>
    dplyr::left_join(n_samples, by = p$food_grouping_variable) |>
    dplyr::mutate(prev = abs_prev / n_all)

  ab_df_sel <- pooled_edges_sample_study |>
    dplyr::mutate(log_ab = log(weight, base = p$log_base)) |>
    dplyr::group_by(
      dplyr::across(dplyr::all_of(c(p$go_down_to, p$food_grouping_variable)))
    ) |>
    dplyr::summarise(
      min_ab       = min(weight),
      max_ab       = max(weight),
      mean_ab      = mean(weight),
      median_ab    = median(weight),
      mad_ab       = mad(weight),
      median_logab = median(log_ab),
      .groups = "drop"
    )

  prev_ab_df_sel <- dplyr::left_join(
    prev_df_sel, ab_df_sel,
    by = c(p$go_down_to, p$food_grouping_variable)
  ) |>
    dplyr::arrange(.data[[p$go_down_to]], .data[[p$food_grouping_variable]])

  list(n_samples = n_samples, prev_ab_df_sel = prev_ab_df_sel)
}
