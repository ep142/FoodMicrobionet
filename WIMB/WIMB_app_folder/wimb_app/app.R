# app.R
# Where Is My Bug? — Shiny dashboard for FoodMicrobionet taxon distribution.
# Requires: bslib, bsicons, shiny, tidyverse, data.table, bipartite, DT, scales
# Place FMBN_plus.rds in the FMBN/ subfolder before launching.

library(shiny)
library(bslib)
library(bsicons)
library(markdown)
library(tidyverse)
library(data.table)
library(DT)
library(scales)

# ---------------------------------------------------------------------------
# Load FMBN once at startup — shared across all sessions
# ---------------------------------------------------------------------------
local_path <- file.path("FMBN", "FMBN_plus.rds")

if (!file.exists(local_path)) {
  stop(
    "FMBN_plus.rds not found in the FMBN/ subfolder.\n",
    "Download it from https://github.com/ep142/FoodMicrobionet and place it ",
    "in the FMBN/ folder next to app.R."
  )
}

FMBN_plus     <- readRDS(local_path)
fmbn_version  <- FMBN_plus$version %||% "unknown"
fmbn_mtime    <- file.mtime(local_path)
fmbn_age_days <- as.numeric(difftime(Sys.time(), fmbn_mtime, units = "days"))
fmbn_stale    <- fmbn_age_days > 365

# ---------------------------------------------------------------------------
# Default parameters (non-UI parameters; UI controls override the others)
# ---------------------------------------------------------------------------
default_p <- list(
  go_down_to             = "species",
  taxon_level            = "genus",
  keep_NA_lowest         = TRUE,
  remove_chloroplast     = FALSE,
  remove_mitochondria    = FALSE,
  remove_no_phylum       = TRUE,
  s_type                 = "Sample",
  food_grouping_variable = "L1",
  min_read_length        = 50,
  run_further_filters    = FALSE,
  gene_targets           = c("ITS1", "ITS2", "ITS1-ITS2",
                              "V1-V3", "V3-V4", "V4"),
  log_base               = 10,
  flip_box_plot          = TRUE,
  violin_min_n           = 10,
  plot_facet_by          = "taxon",
  plot_fill_var          = "rel_prev",
  plot_show_jitter       = TRUE,
  plot_notch             = TRUE,
  dotplot_max_taxa       = 10,
  dotplot_max_foods      = 12,
  dotplot_size_var       = "median_logab",
  min_n_reads            = 5000,
  abbr                   = 15,
  graph_width            = 7,
  graph_height           = 7,
  graph_res              = 150,
  graph_type             = ".jpeg"
)

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- page_sidebar(
  title = "Where Is My Bug?",
  theme = bs_theme(bootswatch = "flatly",
                   base_font  = font_google("Source Sans Pro")),

  tags$head(tags$style(HTML("
    /* --- smaller sidebar text --- */
    .bslib-sidebar-layout .sidebar,
    .bslib-sidebar-layout .sidebar label,
    .bslib-sidebar-layout .sidebar .form-label,
    .bslib-sidebar-layout .sidebar .form-text,
    .bslib-sidebar-layout .sidebar h6 {
      font-size: 0.80rem;
    }
    .bslib-sidebar-layout .sidebar .form-control,
    .bslib-sidebar-layout .sidebar .selectize-input,
    .bslib-sidebar-layout .sidebar .selectize-dropdown {
      font-size: 0.80rem;
    }
    .bslib-sidebar-layout .sidebar .btn {
      font-size: 0.82rem;
    }
  "))),

  sidebar = sidebar(
    width = 380,

    # --- explanatory text ---------------------------------------------------
    accordion(
      open = FALSE,
      accordion_panel(
        "About this app",
        icon = bsicons::bs_icon("info-circle"),
        div(style = "font-size: 0.70rem;", includeMarkdown("explanation.md"))
      )
    ),

    hr(),

    # --- run button (prominent, top of controls) ----------------------------
    actionButton("run", "Run", class = "btn-primary w-100",
                 icon = icon("play")),

    hr(),

    # --- taxonomic controls -------------------------------------------------
    h6("Taxonomic selection", class = "text-muted fw-bold"),

    selectInput(
      "taxon_level", "Search at level",
      choices  = c("genus", "family"),
      selected = "genus"
    ),

    selectInput(
      "go_down_to", "Go down to",
      choices  = c("species", "genus"),
      selected = "species"
    ),

    selectizeInput(
      "taxon_to_search",
      "Taxa (delete default or type to search, max 12)",
      choices  = NULL,
      multiple = TRUE,
      options  = list(
        maxItems    = 12,
        placeholder = "type to search...",
        plugins     = list("remove_button")
      )
    ),

    hr(),

    # --- plot options -------------------------------------------------------
    h6("Plot options", class = "text-muted fw-bold"),

    selectInput(
      "plot_type", "Plot type",
      choices  = c("Boxplot"  = "boxplot",
                   "Violin"   = "violin",
                   "Dot plot" = "dotplot"),
      selected = "boxplot"
    ),

    selectInput(
      "plot_facet_by", "Facet / axis orientation",
      choices  = c("Taxon on facets, food on x" = "food",
                   "Food on facets, taxon on x" = "taxon"),
      selected = "taxon"
    ),

    checkboxInput("plot_show_jitter", "Show jitter points", value = TRUE),
    checkboxInput("plot_notch",       "Notch boxplots",     value = TRUE),
    checkboxInput("flip_box_plot",    "Flip axes",          value = TRUE),

    hr(),

    # --- export options -----------------------------------------------------
    accordion(
      open = FALSE,
      accordion_panel(
        "Export options",
        icon = bsicons::bs_icon("download"),
        textInput("file_name", "File name (no extension)",
                  value = "wimb_plot", placeholder = "e.g. lactococcus_milk"),
        numericInput("graph_width",  "Width (inches)",   value = 7,
                     min = 3, max = 20, step = 0.5),
        numericInput("graph_height", "Height (inches)",  value = 7,
                     min = 3, max = 20, step = 0.5),
        numericInput("graph_res",    "Resolution (dpi)", value = 150,
                     min = 72, max = 600, step = 50),
        selectInput("graph_type", "Format",
                    choices  = c("JPEG" = ".jpeg", "PNG" = ".png",
                                 "PDF"  = ".pdf",  "SVG" = ".svg"),
                    selected = ".jpeg")
      )
    ),

    hr(),

    # --- FMBN version badge -------------------------------------------------
    tags$small(
      class = "text-muted",
      tags$b("FoodMicrobionet "), fmbn_version,
      tags$br(),
      sprintf("Local cache: %.0f days old", fmbn_age_days),
      if (fmbn_stale) tags$span(
        " \u26a0 consider updating",
        style = "color: orange; font-weight: bold;"
      )
    )
  ),

  # --- run log (compact card, always above the plot tabs) ------------------
  card(
    fill = FALSE,
    card_header(
      class = "py-1 d-flex align-items-center gap-2",
      style = "font-size: 0.78rem;",
      bsicons::bs_icon("terminal"), "Run log"
    ),
    card_body(
      class = "p-2",
      div(
        style = "height: 90px; overflow-y: auto;",
        uiOutput("log_pane")
      )
    )
  ),

  # --- main panel -----------------------------------------------------------
  navset_card_tab(

    nav_panel(
      "Distribution plot",
      icon = bsicons::bs_icon("bar-chart"),
      card(
        card_body(
          uiOutput("run_status"),
          plotOutput("distplot", height = "550px")
        ),
        card_footer(
          downloadButton("dl_distplot", "Download plot",
                         class = "btn-sm btn-outline-secondary")
        )
      )
    ),

    nav_panel(
      "Dot plot",
      icon = bsicons::bs_icon("grid-3x3-gap"),
      card(
        card_body(
          plotOutput("dotplot", height = "550px")
        ),
        card_footer(
          downloadButton("dl_dotplot", "Download plot",
                         class = "btn-sm btn-outline-secondary")
        )
      )
    ),

    nav_panel(
      "Summary table",
      icon = bsicons::bs_icon("table"),
      card(
        card_body(
          DT::dataTableOutput("prev_ab_table")
        ),
        card_footer(
          downloadButton("dl_prev_ab", "Download (TSV)",
                         class = "btn-sm btn-outline-secondary")
        )
      )
    ),

    nav_panel(
      "Abbreviations",
      icon = bsicons::bs_icon("card-text"),
      card(
        card_body(
          DT::dataTableOutput("abbr_table")
        ),
        card_footer(
          downloadButton("dl_abbr", "Download (TSV)",
                         class = "btn-sm btn-outline-secondary")
        )
      )
    ),

    nav_panel(
      "Bibliography",
      icon = bsicons::bs_icon("journal-text"),
      card(
        card_body(
          DT::dataTableOutput("bib_table")
        ),
        card_footer(
          downloadButton("dl_bib", "Download (TSV)",
                         class = "btn-sm btn-outline-secondary")
        )
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  # --- persistent message store ---------------------------------------------
  app_messages <- reactiveVal(character(0))

  # --- update go_down_to choices based on taxon_level ----------------------
  observeEvent(input$taxon_level, {
    choices <- if (input$taxon_level == "family") {
      c("genus")
    } else {
      c("species", "genus")
    }
    updateSelectInput(session, "go_down_to",
                      choices = choices, selected = choices[1])
  })

  # --- server-side selectize: update taxon list reactively -----------------
  observeEvent(input$taxon_level, {
    choices <- FMBN_plus$taxa |>
      dplyr::filter(!is.na(.data[[input$taxon_level]])) |>
      dplyr::pull(.data[[input$taxon_level]]) |>
      unique() |>
      sort()
    # pre-select a sensible default so the field is never blank on startup
    default_taxon <- if (input$taxon_level == "genus") "Lactobacillus"
                     else                               "Lactobacillaceae"
    selected <- if (default_taxon %in% choices) default_taxon else character(0)
    updateSelectizeInput(session, "taxon_to_search",
                         choices  = choices,
                         server   = TRUE,
                         selected = selected)
  }, ignoreNULL = FALSE)

  # --- reactive data: only fires on Run button ------------------------------
  results <- eventReactive(input$run, {
    req(input$taxon_to_search)

    if (length(input$taxon_to_search) > 12) {
      showNotification("Please select 12 or fewer taxa.", type = "warning")
      return(NULL)
    }

    # start collecting messages for this run (log accumulates across runs)
    collected_msgs <- character(0)

    # assemble p from UI inputs and defaults
    p <- modifyList(default_p, list(
      taxon_level      = input$taxon_level,
      go_down_to       = input$go_down_to,
      taxon_to_search  = input$taxon_to_search,
      plot_type        = input$plot_type,
      plot_facet_by    = input$plot_facet_by,
      plot_show_jitter = input$plot_show_jitter,
      plot_notch       = input$plot_notch,
      flip_box_plot    = input$flip_box_plot,
      graph_width      = input$graph_width,
      graph_height     = input$graph_height,
      graph_res        = input$graph_res,
      graph_type       = input$graph_type
    ))

    result <- withCallingHandlers(
      withProgress(message = "Loading data...", value = 0, {

        setProgress(0.2, detail = "Filtering taxa and edges")
        edges_result <- tryCatch(
          build_pooled_edges(FMBN_plus, p),
          error = function(e) {
            showNotification(conditionMessage(e), type = "error",
                             duration = 15)
            NULL
          }
        )
        if (is.null(edges_result)) return(NULL)

        setProgress(0.6, detail = "Computing prevalence and abundance")
        pa_result <- build_prev_ab(
          edges_result$pooled_edges_sample_study,
          edges_result$filt_samples,
          p
        )

        setProgress(0.9, detail = "Building plot data")

        pes <- edges_result$pooled_edges_sample_study
        df_to_plot <- pes |>
          dplyr::mutate(
            taxon_abbr      = stringr::str_trunc(
              .data[[p$go_down_to]],
              width = p$abbr, side = "center"),
            food_group_abbr = stringr::str_trunc(
              .data[[p$food_grouping_variable]],
              width = p$abbr, side = "center"),
            log_ab          = log(weight, base = p$log_base)
          )

        species_levels   <- sort(unique(pes[[p$go_down_to]]))
        taxa_abbr_levels <- sort(unique(df_to_plot$taxon_abbr))

        df_to_plot <- df_to_plot |>
          dplyr::mutate(
            taxon_abbr = factor(taxon_abbr, levels = taxa_abbr_levels),
            !!rlang::sym(p$go_down_to) := factor(
              .data[[p$go_down_to]], levels = species_levels)
          ) |>
          dplyr::left_join(
            dplyr::select(
              pa_result$prev_ab_df_sel,
              dplyr::all_of(c(p$go_down_to, p$food_grouping_variable)),
              rel_prev = prev
            ),
            by = c(p$go_down_to, p$food_grouping_variable)
          )

        setProgress(1)

        list(
          p                         = p,
          df_to_plot                = df_to_plot,
          prev_ab_df_sel            = pa_result$prev_ab_df_sel,
          pooled_edges_sample_study = edges_result$pooled_edges_sample_study,
          filt_studies              = edges_result$filt_studies,
          taxon_label               = edges_result$taxon_label,
          nstudies_wtarget          = edges_result$nstudies_wtarget,
          nsamples_wtarget          = edges_result$nsamples_wtarget,
          nstudies_total            = edges_result$nstudies_total,
          nsamples_total            = edges_result$nsamples_total
        )
      }),

      # intercept all message() and warning() calls anywhere in the pipeline
      message = function(m) {
        collected_msgs <<- c(collected_msgs,
                             stringr::str_trim(conditionMessage(m)))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        collected_msgs <<- c(collected_msgs,
                             paste("\u26a0",
                                   stringr::str_trim(conditionMessage(w))))
        invokeRestart("muffleWarning")
      }
    )

    # append this run's messages to the permanent log, with a timestamp header
    run_header <- paste0("── ", format(Sys.time(), "%H:%M:%S"), " · ",
                         paste(input$taxon_to_search, collapse = ", "), " ──")
    app_messages(c(isolate(app_messages()), run_header, collected_msgs))
    result
  })

  # --- run log pane ---------------------------------------------------------
  output$log_pane <- renderUI({
    msgs <- app_messages()
    if (length(msgs) == 0) {
      return(tags$span(class = "text-muted",
                       style = "font-size:0.78rem;",
                       "No messages yet — click Run to start."))
    }
    tags$pre(
      class  = "mb-0",
      style  = "font-size:0.75rem; white-space:pre-wrap; word-break:break-word;
                background:transparent; border:none; padding:0;",
      paste(msgs, collapse = "\n")
    )
  })

  # --- status message -------------------------------------------------------
  output$run_status <- renderUI({
    r <- results()
    req(r)
    tags$p(
      class = "text-muted small mb-1",
      sprintf(
        "%s found in %d studies and %d samples (out of %d studies, %d samples in FMBN).",
        r$taxon_label,
        r$nstudies_wtarget, r$nsamples_wtarget,
        r$nstudies_total,   r$nsamples_total
      )
    )
  })

  # --- abbreviations table --------------------------------------------------
  abbr_data <- reactive({
    r <- results()
    req(r)
    taxon_res <- resolve_abbr(
      sort(unique(as.character(r$df_to_plot[[r$p$go_down_to]]))),
      r$p$abbr
    )
    food_res <- resolve_abbr(
      sort(unique(r$df_to_plot[[r$p$food_grouping_variable]])),
      r$p$abbr
    )
    rbind(
      data.frame(
        Type          = r$p$go_down_to,
        `Full label`  = sort(unique(as.character(r$df_to_plot[[r$p$go_down_to]]))),
        `Plot label`  = taxon_res$abbr,
        check.names   = FALSE
      ),
      data.frame(
        Type          = r$p$food_grouping_variable,
        `Full label`  = sort(unique(r$df_to_plot[[r$p$food_grouping_variable]])),
        `Plot label`  = food_res$abbr,
        check.names   = FALSE
      )
    )
  })

  output$abbr_table <- DT::renderDataTable({
    DT::datatable(
      abbr_data(),
      rownames = FALSE,
      filter   = "top",
      options  = list(pageLength = 30, scrollX = TRUE),
      caption  = "Mapping between full labels and abbreviated plot labels"
    )
  })

  output$dl_abbr <- downloadHandler(
    filename = function() paste0(input$file_name, "_abbr.tsv"),
    content  = function(file) readr::write_tsv(abbr_data(), file)
  )

  # --- weightlabel helper ---------------------------------------------------
  weightlabel <- reactive({
    r <- results()
    req(r)
    if (r$p$log_base == 2) "log2(rel. abundance)" else "log10(rel. abundance)"
  })

  # --- distribution plot ----------------------------------------------------
  # Returns list(plt, msgs) so we can propagate messages without a circular
  # reactive dependency (reading app_messages inside a reactive that writes it).
  distplot_obj <- reactive({
    r <- results()
    req(r)
    if (r$p$plot_type == "dotplot") return(list(plt = NULL, msgs = character(0)))

    new_msgs <- character(0)
    plt <- withCallingHandlers(
      make_distplot(r$df_to_plot, r$p, weightlabel()),
      message = function(m) {
        new_msgs <<- c(new_msgs, stringr::str_trim(conditionMessage(m)))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        new_msgs <<- c(new_msgs,
                       paste("\u26a0", stringr::str_trim(conditionMessage(w))))
        invokeRestart("muffleWarning")
      }
    )
    list(plt = plt, msgs = new_msgs)
  })

  # Append any plot-level messages to the persistent store.
  # Using observeEvent (not reactive) avoids a circular dependency.
  observeEvent(distplot_obj(), {
    obj <- distplot_obj()
    if (!is.null(obj) && length(obj$msgs) > 0) {
      app_messages(c(isolate(app_messages()), obj$msgs))
    }
  }, ignoreNULL = TRUE)

  output$distplot <- renderPlot({
    r <- results()
    req(r)
    if (r$p$plot_type == "dotplot") {
      make_dotplot(r$prev_ab_df_sel, r$p, weightlabel())
    } else {
      distplot_obj()$plt
    }
  })

  output$dl_distplot <- downloadHandler(
    filename = function() {
      r <- results()
      ext <- if (!is.null(r)) r$p$graph_type else ".jpeg"
      paste0(input$file_name, ext)
    },
    content = function(file) {
      r <- results()
      req(r)
      plt <- if (r$p$plot_type == "dotplot") {
        make_dotplot(r$prev_ab_df_sel, r$p, weightlabel())
      } else {
        distplot_obj()$plt
      }
      ggplot2::ggsave(file, plt,
                      width  = r$p$graph_width,
                      height = r$p$graph_height,
                      dpi    = r$p$graph_res)
    }
  )

  # --- dot plot -------------------------------------------------------------
  dotplot_obj <- reactive({
    r <- results()
    req(r)
    make_dotplot(r$prev_ab_df_sel, r$p, weightlabel())
  })

  output$dotplot <- renderPlot({
    dotplot_obj()
  })

  output$dl_dotplot <- downloadHandler(
    filename = function() {
      r <- results()
      ext <- if (!is.null(r)) r$p$graph_type else ".jpeg"
      paste0(input$file_name, "_dot", ext)
    },
    content = function(file) {
      r <- results()
      req(r)
      ggplot2::ggsave(file, dotplot_obj(),
                      width  = r$p$graph_width,
                      height = r$p$graph_height,
                      dpi    = r$p$graph_res)
    }
  )

  # --- summary table --------------------------------------------------------
  output$prev_ab_table <- DT::renderDataTable({
    r <- results()
    req(r)
    DT::datatable(
      r$prev_ab_df_sel,
      rownames = FALSE,
      filter   = "top",
      options  = list(pageLength = 15, scrollX = TRUE),
      caption  = paste("Prevalence and abundance summary for", r$taxon_label)
    ) |>
      DT::formatRound(
        columns = c("prev", "min_ab", "max_ab", "mean_ab",
                    "median_ab", "mad_ab", "median_logab"),
        digits  = 3
      )
  })

  output$dl_prev_ab <- downloadHandler(
    filename = function() paste0(input$file_name, "_summary.tsv"),
    content  = function(file) {
      r <- results()
      req(r)
      readr::write_tsv(r$prev_ab_df_sel, file)
    }
  )

  # --- bibliography table ---------------------------------------------------
  bib_data <- reactive({
    r <- results()
    req(r)
    my_studies <- dplyr::distinct(r$pooled_edges_sample_study, studyId) |>
      dplyr::pull(studyId)
    r$filt_studies |>
      dplyr::filter(studyId %in% my_studies) |>
      dplyr::select(studyId, target, region, bioproject,
                    food_group, DOI, ref_complete)
  })

  output$bib_table <- DT::renderDataTable({
    DT::datatable(
      bib_data(),
      rownames = FALSE,
      filter   = "top",
      escape   = FALSE,
      options  = list(pageLength = 15, scrollX = TRUE),
      caption  = "FoodMicrobionet studies contributing data for selected taxa"
    )
  })

  output$dl_bib <- downloadHandler(
    filename = function() paste0(input$file_name, "_bib.tsv"),
    content  = function(file) {
      readr::write_tsv(bib_data(), file)
    }
  )
}

# ---------------------------------------------------------------------------
shinyApp(ui, server)
