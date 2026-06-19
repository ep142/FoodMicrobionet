# make_dotplot.R
# Dot plot function for the Where Is My Bug? app.
# Sourced automatically by Shiny before app.R runs.

make_dotplot <- function(df_summary, p, weightlabel) {

  taxon_var <- p$go_down_to
  food_var  <- p$food_grouping_variable

  # enforce caps
  top_taxa <- df_summary |>
    dplyr::group_by(.data[[taxon_var]]) |>
    dplyr::summarise(overall_med = median(median_logab, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::slice_max(overall_med, n = p$dotplot_max_taxa,
                     with_ties = FALSE) |>
    dplyr::pull(.data[[taxon_var]])

  top_foods <- df_summary |>
    dplyr::group_by(.data[[food_var]]) |>
    dplyr::summarise(n_taxa = dplyr::n_distinct(.data[[taxon_var]]),
                     .groups = "drop") |>
    dplyr::slice_max(n_taxa, n = p$dotplot_max_foods,
                     with_ties = FALSE) |>
    dplyr::pull(.data[[food_var]])

  df_plot <- df_summary |>
    dplyr::filter(.data[[taxon_var]] %in% top_taxa,
                  .data[[food_var]]  %in% top_foods)

  if (nrow(df_plot) == 0) {
    warning("make_dotplot: no data after applying caps.")
    return(invisible(NULL))
  }

  # abbreviations using p$abbr — widen if needed to avoid label collisions
  taxon_abbr_result <- resolve_abbr(unique(df_plot[[taxon_var]]), p$abbr)
  food_abbr_result  <- resolve_abbr(unique(df_plot[[food_var]]),  p$abbr)

  if (taxon_abbr_result$width != p$abbr)
    message(sprintf("Dot plot: taxon label width increased from %d to %d to avoid conflicts.",
                    p$abbr, taxon_abbr_result$width))
  if (food_abbr_result$width != p$abbr)
    message(sprintf("Dot plot: food label width increased from %d to %d to avoid conflicts.",
                    p$abbr, food_abbr_result$width))

  taxon_map <- stats::setNames(taxon_abbr_result$abbr, unique(df_plot[[taxon_var]]))
  food_map  <- stats::setNames(food_abbr_result$abbr,  unique(df_plot[[food_var]]))

  df_plot <- df_plot |>
    dplyr::mutate(
      taxon_abbr      = taxon_map[.data[[taxon_var]]],
      food_group_abbr = food_map[.data[[food_var]]]
    )

  # axis assignment
  if (p$plot_facet_by == "food") {
    x_var   <- "taxon_abbr"
    y_var   <- "food_group_abbr"
    x_label <- "taxon"
    y_label <- p$food_grouping_variable
  } else {
    x_var   <- "food_group_abbr"
    y_var   <- "taxon_abbr"
    x_label <- p$food_grouping_variable
    y_label <- "taxon"
  }

  # size variable and scale
  size_var   <- p$dotplot_size_var
  size_label <- if (size_var == "median_logab") weightlabel else "median rel. abundance"

  if (size_var == "median_logab") {
    size_scale <- ggplot2::scale_size_continuous(name = size_label, range = c(2, 10))
  } else {
    size_scale <- ggplot2::scale_size_continuous(
      name   = size_label,
      range  = c(2, 10),
      trans  = scales::log_trans(base = p$log_base),
      labels = scales::label_scientific()
    )
  }

  # build plot
  if (p$plot_fill_var == "rel_prev") {
    plt <- ggplot2::ggplot(df_plot,
                           ggplot2::aes(x = .data[[x_var]],
                                        y = .data[[y_var]],
                                        size = .data[[size_var]])) +
      ggplot2::geom_point(ggplot2::aes(fill = prev),
                          shape = 21, alpha = 0.85) +
      ggplot2::scale_fill_viridis_c(option = "plasma", direction = -1,
                                    name = "rel. prev.")
  } else if (p$plot_fill_var == "spoilage") {
    warning("make_dotplot: spoilage fill not available for dot plot. ",
            "Falling back to plain fill.", call. = FALSE)
    plt <- ggplot2::ggplot(df_plot,
                           ggplot2::aes(x = .data[[x_var]],
                                        y = .data[[y_var]],
                                        size = .data[[size_var]])) +
      ggplot2::geom_point(shape = 21, fill = "steelblue", alpha = 0.85)
  } else {
    plt <- ggplot2::ggplot(df_plot,
                           ggplot2::aes(x = .data[[x_var]],
                                        y = .data[[y_var]],
                                        size = .data[[size_var]])) +
      ggplot2::geom_point(shape = 21, fill = "steelblue", alpha = 0.85)
  }

  plt <- plt +
    size_scale +
    ggplot2::labs(x = x_label, y = y_label) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 90, hjust = 1,
                                               vjust = 0.5),
      panel.grid.major = ggplot2::element_line(colour = "grey90"),
      legend.position  = "right"
    )

  if (p$flip_box_plot) plt <- plt + ggplot2::coord_flip()

  plt
}
