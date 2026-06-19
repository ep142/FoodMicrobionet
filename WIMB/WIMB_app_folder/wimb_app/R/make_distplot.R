# make_distplot.R
# Box / violin plot function for the Where Is My Bug? app.
# Sourced automatically by Shiny before app.R runs.

fermspoil_scale <- c("blue", "red", "cyan", "mediumorchid")
names(fermspoil_scale) <- c("Unspoiled", "Spoiled", "Fermented", "Fermented+Spoiled")

make_distplot <- function(df, p, weightlabel) {

  if (p$plot_facet_by == "food") {
    x_var          <- "taxon_abbr"
    x_full_col     <- p$go_down_to
    facet_var      <- "food_group_abbr"
    facet_full_col <- p$food_grouping_variable
    x_label        <- "taxon"
  } else {
    x_var          <- "food_group_abbr"
    x_full_col     <- p$food_grouping_variable
    facet_var      <- "taxon_abbr"
    facet_full_col <- p$go_down_to
    x_label        <- p$food_grouping_variable
  }

  # Guard against abbreviation collisions in the facet variable.
  # E.g. "Lactococcus hircilactis", "Lactococcus lactis", and
  # "Lactococcus raffinolactis" all truncate to "Lactoc...lactis" at width 15,
  # which merges their rows and breaks the fill colour mapping.
  # Widen the abbreviation one character at a time until labels are unique.
  facet_abbr_res <- resolve_abbr(unique(df[[facet_full_col]]), p$abbr)
  if (facet_abbr_res$width != p$abbr) {
    message(sprintf(
      "Label width increased from %d to %d to avoid abbreviation conflicts in facet labels.",
      p$abbr, facet_abbr_res$width))
    df <- df |>
      dplyr::mutate(!!facet_var := stringr::str_trunc(
        .data[[facet_full_col]], facet_abbr_res$width, side = "center"))
  }

  # single-taxon guard: drop facet when only one level in the data
  n_facets <- dplyr::n_distinct(df[[facet_var]])

  # --- warn when relative prevalence is missing for some groups -----------
  n_na_prev <- sum(is.na(df$rel_prev))
  if (n_na_prev > 0) {
    warning(
      sprintf("%d observation(s) have no relative prevalence (shown in light orange). ",
              n_na_prev),
      "This usually means that food-group totals could not be computed for ",
      "those taxon \u00d7 food-group combinations.",
      call. = FALSE
    )
  }

  # --- decide effective plot type ------------------------------------------
  effective_type <- p$plot_type
  if (p$plot_type == "violin") {
    group_sizes <- df |>
      dplyr::group_by(.data[[x_var]], .data[[facet_var]]) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
      dplyr::pull(n)
    
    # fall back only if the MAJORITY of groups are too small,
    # not just because one sparse combination exists
    pct_too_small <- mean(group_sizes < p$violin_min_n)
    if (pct_too_small > 0.5) {
      effective_type <- "boxplot"
      message("make_distplot: more than half of groups have fewer than ",
              p$violin_min_n, " observations — falling back to boxplot.")
    }
  }

  outlier_shape <- if (p$plot_show_jitter) NA else 19

  plt <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]], y = log_ab))

  # jitter layer
  if (p$plot_show_jitter) {
    if (p$plot_fill_var == "spoilage") {
      plt <- plt +
        ggplot2::geom_jitter(ggplot2::aes(colour = spoilage),
                             alpha = 0.4, size = 1, width = 0.1, shape = 16) +
        ggplot2::scale_colour_manual(values = fermspoil_scale,
                                     na.value = "grey80")
    } else {
      plt <- plt +
        ggplot2::geom_jitter(colour = "grey50", alpha = 0.25,
                             size = 0.8, width = 0.15)
    }
  }

  # main geom
  if (effective_type == "violin") {
    if (p$plot_fill_var == "rel_prev") {
      plt <- plt +
        ggplot2::geom_violin(ggplot2::aes(fill = rel_prev), alpha = 0.7,
                             draw_quantiles = c(0.25, 0.5, 0.75),
                             scale = "width") +
        ggplot2::scale_fill_viridis_c(option = "plasma", direction = -1,
                                      name = "rel. prev.",
                                      na.value = "#FFCC99")
    } else {
      plt <- plt +
        ggplot2::geom_violin(alpha = 0.7,
                             draw_quantiles = c(0.25, 0.5, 0.75),
                             scale = "width")
    }
  } else {
    if (p$plot_fill_var == "rel_prev") {
      plt <- plt +
        ggplot2::geom_boxplot(ggplot2::aes(fill = rel_prev), alpha = 0.7,
                              notch = p$plot_notch,
                              outlier.shape = outlier_shape) +
        ggplot2::scale_fill_viridis_c(option = "plasma", direction = -1,
                                      name = "rel. prev.",
                                      na.value = "#FFCC99")
    } else {
      plt <- plt +
        ggplot2::geom_boxplot(alpha = 0.7, notch = p$plot_notch,
                              outlier.shape = outlier_shape)
    }
  }

  # facet only when more than one level
  if (n_facets > 1) {
    plt <- plt + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)))
  }

  plt <- plt +
    ggplot2::labs(x = x_label, y = weightlabel) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)
    )

  if (p$log_base == 10) {
    plt <- plt + ggplot2::annotation_logticks(sides = "l", alpha = 0.5)
  }
  if (p$flip_box_plot) plt <- plt + ggplot2::coord_flip()

  plt
}
