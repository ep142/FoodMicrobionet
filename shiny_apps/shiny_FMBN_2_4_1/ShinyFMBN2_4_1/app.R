# ShinyFMBN v2_4_1-----------------------------------------------------

# a Shiny app to explore, filter and extract data from FoodMicrobionet (v4.1, compatible with 3.2.8)

# preparatory steps -------------------------------------------------------

# install/load packages ---------------------------------------------------

.cran_packages <-
  c(
    "shiny",
    "DT",
    "tidyverse",
    "reshape2",
    "stringr",
    "lubridate",
    "igraph",
    "magrittr",
    "randomcoloR",
    "forcats"
  )
.bioc_packages <- c(
                    "BiocManager", 
                    "phyloseq"
                    )
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  if(!.inst[1]) {
    install.packages("BiocManager")
    .inst <- .bioc_packages %in% installed.packages()
  }
  if(any(!.inst[2:length(.inst)])) {
    BiocManager::install(.bioc_packages[!.inst], ask = F)
  }
}
sapply(c(.cran_packages, .bioc_packages), require,
       character.only = TRUE)

# options -----------------------------------------------------------------

# debug mode
# set to true and extra information for cross checking will appear in the filter and aggregate tabs
debug_mode <- F

# options for filtering
# set to 0 so the user has to move them
minprev <-
  0 # the default minimum for prevalence filtering (see Export pane/Prevalence graph tab)
minab <-
  0 # the default minimum for abundance filtering (see Export pane/Prevalence graph tab)

# default options for ggsave (changes here will apply to every ggsave in the app)
g_width <- 7
g_height <- 5
g_units = "in"
g_dpi <- 150
g_ext <- "jpg" # use "jpg", "tif" or "pdf"

# max taxa for barplots and boxplots
max_taxa <- 25
# max "samples" or sample categories for barplots and boxplots
max_samples <- 25

# nchar for str_wrap in graphs, needed to fine tune x-axis labels on some graphs
str_wrap_length <- 20

# renumber samples in lists for export
renumber_samples <- FALSE


# the random palette
rpalette <- distinctColorPalette(max_taxa)


# load and assemble  -------------------------------------
# load the list (from the data folder) using a platform independent path

FMBN <- readRDS(file.path("data", "FMBN.rds"))

studies <- FMBN$studies
samples <- FMBN$samples
edges <- FMBN$edges
taxa <- FMBN$taxa
references <- FMBN$references
version <- FMBN$version


# generate minitables for display
studies_disp <- studies %>%
  dplyr::select(
    studyId,
    FMBN_version,
    target,
    region,
    platform,
    seq_center,
    bioinf_software,
    OTU_picking,
    assign_tax_method,
    tax_database,
    Seq_accn,
    Seq_accn_link,
    short_descr,
    samples,
    ref_short,
    DOI,
    DOI_link
  )
studies_disp <- studies_disp %>%
  mutate(
    target = paste(target, ", ", region, sep = ""),
    pipeline = str_c(
      bioinf_software,
      OTU_picking,
      assign_tax_method,
      tax_database,
      sep = "; "
    ),
    SRA_ENA = ifelse(
      is.na(Seq_accn),
      "not available",
      paste0(
        "<a href='",
        Seq_accn_link,
        "' target='_blank'>",
        Seq_accn,
        "</a>"
      )
    ),
    DOIlink = ifelse(
      DOI == "unpublished data",
      "unpublished data",
      paste0("<a href='",
             DOI_link,
             "' target='_blank'>",
             DOI, "</a>")
    )
  ) %>%
  dplyr::select(
    studyId,
    FMBN_version,
    target,
    platform,
    pipeline,
    samples,
    SRA_ENA,
    ref_short,
    DOIlink,
    short_descr
  )

# generate lists for check box groups
s_type_list <- unique(samples$s_type)
nature_list <- unique(samples$nature)
process_list <- unique(samples$process)
spoilage_list <- unique(samples$spoilage)
target1_list <- c(DNA = "16S_DNA", RNA = "16S_RNA")
target2_list <- unique(samples$target2)

# generates a transformed table for sample selection
ssamples <- samples %>%
  mutate_at(c("s_type", "nature", "process", "spoilage", "target1", "target2"),
            as.factor)

# generate (non reactive) summary data on studies, samples, taxa ---------------

n_studies <- nrow(studies)
n_samples <- nrow(samples)
n_food_groups <- n_distinct(samples$L1)
n_foodId <- n_distinct(samples$foodId)
n_llabel <- n_distinct(samples$llabel)
n_taxa <- nrow(taxa)

FMBN_summary_text <-
  paste(
    "There are",
    n_studies,
    "studies and",
    n_samples,
    "samples in this version of FoodMicrobionet.",
    "The samples belong to",
    n_food_groups,
    "major food groups and",
    n_foodId,
    "different foods.",
    "There are",
    n_llabel,
    "different combinations of food, nature, process, fermentation/spoilage.",
    n_taxa,
    "taxa have been identified at different taxonomic levels."
  )


# user interface ----------------------------------------------------------

ui <- navbarPage(
  "Shiny FMBN",
  
  # View tab --------------------------------------------------------------
  
  tabPanel("Explore",
           fluidPage(
             # row for studies table
             fluidRow(
               column(3,
                      h4("FMBN studies"),
                      p("Explore studies in FMBN by navigating the table on the right."),
                      tags$div(tags$ul(
                        tags$li(
                          tags$span("Use the search box to find studies by keywords in any field")
                        ),
                        tags$li(
                          tags$span("Use the links to access studies on SRA/ENA or articles via DOI.")
                        ),
                        tags$li(
                          tags$span(
                            "Hover with the mouse over a pipeline or description cell to get a full description."
                          )
                        )
                      )
                      )
               ),
               column(9, 
                      div(DT::dataTableOutput("studies_table"), style = "font-size:70%")
               )
             ),
             # end row for studies table
             
             # row for summary stats on FMBN
             fluidRow(column(3,
                             hr(),
                             h4("Summary statistics on FMBN")
             ),
             column(9,
                    hr(),
                    p(paste(version, FMBN_summary_text))
             )
             ),
             # end row for summary stats
             
             # row for selecting and viewing a study
             fluidRow(
               column(3,
                      hr(),
                      h4("Samples"),
                      uiOutput('choose_study')
               ),
               column(9,
                      hr(),
                      div(DT::dataTableOutput("vsamples_table"), style = "font-size:70%")
               )
             )
             # end row sol selecting and viewing a study
             
           )
           # end fluid page View tab
           
  ),
  # end View tab
  
  # Filter tab --------------------------------------------------
  tabPanel("Filter",
           fluidPage(
             
             # Filter instructions  and widgets ----------------------------
             fluidRow(
               column(2,
                      selectInput(
                        "select_study",
                        label = "Select a subset of studies",
                        choices = studies$studyId,
                        selected = NULL,
                        multiple = T
                      ),
                      checkboxInput("study_ch_box", label = "apply", value = F),
                      uiOutput('choose_food_group'),
                      checkboxInput("food_group_ch_box", label = "apply", value = F),
                      uiOutput('choose_foodId'),
                      checkboxInput("foodId_ch_box", label = "apply", value = F),
                      hr(),
                      strong("Instructions"),
                      p(
                        "Use dropdown menus to select by study, food group, food code (FoodEx2). Selected samples will appear in the table, which can be used as a guidance for refining your selection.",
                        "Dropdown filters are applied in the order they appear in this page. Use the <apply> checkbox to apply/remove a filter from a drop down menu. ",
                        "When you are done, you can use the filters at the bottom of the table (sample type, nature, process, target) to further refine your search.",
                        style = "font-size:80%"
                      )
               ),
               # end left column Filter tab
               
               # Filtered table ----------------------------------------------
               column(10,
                      div(
                        DT::dataTableOutput("ssamplestable"), style = "font-size:80%"
                      )
               )
             ),
             # end row Filter pane
             
             # diagnostic for the filter pane -----------------------------
             # will appear if debug_mode == T (line 92)
             fluidRow(column(2,
                             if (debug_mode) {
                               h5("Sample ids, table in this tab")
                             }),
                      column(10,
                             if (debug_mode) {
                               h5("Sample labels, table in this tab")
                             })),
             fluidRow(column(2,
                             if (debug_mode) {
                               textOutput("samples_ids")
                             }),
                      column(10,
                             if (debug_mode) {
                               textOutput("samples_label")
                             }))
           )
  ),
  # end Filter tab
  
  # Aggregate tab -----------------------------------------------------------
  tabPanel("Aggregate",
           fluidPage(
             fluidRow(
               column(2,
                      selectInput(
                        "sample_agg_level",
                        label = "Select the sample aggregation level",
                        choices = c("sample", "exp. code"),
                        selected = "sample",
                        multiple = F
                      ),
                      selectInput(
                        "tax_agg_level",
                        label = "Select the taxonomic aggregation level",
                        choices = c("species", "genus", "family", "class"),
                        selected = "species",
                        multiple = F
                      ),
                      selectInput(
                        "show_what",
                        label = "Show",
                        choices = c("original", "aggregated"),
                        selected = "original",
                        multiple = F
                      ),
                      hr(),
                      strong("Instructions"),
                      p(
                        "Once you have selected a subset of samples in the Filter tab, ",
                        "tables for studies, samples, taxa (in descending order of ",
                        "relative abundance), and references will appear on the left.",
                        style = "font-size:80%"
                      ),
                      p(
                        "If you are happy with your selection, use the tools in this tab ",
                        "to (optionally) aggregate samples and taxa and to show original or ",
                        "aggregated data and then move to the Export tab to save your files. ",
                        "Otherwise go back to the Filter tab and start over.",
                        style = "font-size:80%"
                      ),
                      strong("Please note:"),
                      p(
                        "When combining studies using different platforms/pipelines, you should choose",
                        "genus or above as taxonomic aggregation level.",
                        style = "font-size:80%"
                      ),
                      p(
                        "When you choose the expanded food code as a sample aggregation level, ",
                        "consider the risks related to aggregating samples for which different",
                        "wet and dry lab procedure and bionformatics pipelines were used. ",
                        "In addition at the latter aggregation level a heuristic is used to ",
                        "recalculate the number of sequences for each group.",
                        style = "font-size:80%"
                      ),
                      p(
                        "Refreshing the display for large tables may take some time.",
                        style = "font-size:80%"
                      )
               ),
               # end left column Aggregate tab
               
               column(10,
                      textOutput("subset_info"),
                      h4("Your studies"),
                      div(DT::dataTableOutput("f_studies_table"), style = "font-size:80%"),
                      h4("Your samples"),
                      div(DT::dataTableOutput("f_samples_table"), style = "font-size:80%"),
                      h4("Your taxa"),
                      div(DT::dataTableOutput("f_taxa_table"), style = "font-size:80%"),
                      h4("Your references"),
                      div(DT::dataTableOutput("f_refs"), style = "font-size:80%")
               )
               #end right column Aggregate tab
             ),
             
             fluidRow(column(2,
                             # diagnostic, will appear if debug_mode == TRUE at line 92
                             if (debug_mode) {
                               h5("row numbers from filtered table") # it is here for diagnostic purposes
                             }),
                      column(10,
                             # diagnostic, will appear if debug_mode == TRUE at line 92)
                             if (debug_mode) {
                               h5("sliced table from samples df")
                             })),
             fluidRow(column(2,
                             # diagnostic, will appear if debug_mode == TRUE at line 92
                             if (debug_mode) {
                               textOutput("sel_row_numbers") # it is here for diagnostic purposes
                             },
                             if (debug_mode) {
                               textOutput("no_agg") # it is here for diagnostic purposes
                             }
                             ),
                      column(10,
                             # diagnostic, will appear if debug_mode == TRUE at line 92)
                             if (debug_mode) {
                               DT::dataTableOutput("agg_samples_mt")
                             })
             )
           )
           # end fluidPage Aggregate tab
  ),
  # end Aggregate tab
  
  # Export tab --------------------------------------------------------------
  tabPanel("Export",
           sidebarLayout(
             # Sidebar panel for inputs ----
             sidebarPanel(
               h4("Show your summary"),
               p(
                 "Calculating outputs may take time for large datasets and may freeze the display. ",
                 "You have to select the <show me!> checkbox to show your summaries in the tabs of the main panel.",
                 style = "font-size:80%"
               ),
               checkboxInput("show_summary", "show me!", value = F),
               if (debug_mode)
                 verbatimTextOutput("show_output"),
               hr(),
               h4("Export your files."),
               p(
                 "Use the text input to provide a prefix which will be included ",
                 "in all your filenames.",
                 style = "font-size:80%"
               ),
               textInput("fn_prefix", "Filename prefix", value = "myFMBN"),
               p(
                 "Use the action buttons to export files extracted from FoodMicrobionet in the output folder:",
                 style = "font-size:80%"
               ),
               actionButton("miniFMBN", "all"),
               actionButton("aggFMBN", "agg"),
               actionButton("phyloseq", "phy"),
               actionButton("conet", "con"),
               actionButton("gml", "gml"),
               actionButton("ref", "ref"),
               hr(),
               tags$div(
                 tags$ul(
                   tags$li(
                     tags$span(
                       "all: a .rds file with filtered FMBN tables, no aggregation (subfolder minifmbn)"
                     )
                   ),
                   tags$li(
                     tags$span(
                       "agg: a .rds file with FMBN tables after aggregation, for further graphical and statistical analysis (subfolder aggdata)"
                     )
                   ),
                   tags$li(
                     tags$span(
                       "phy: a phyloseq class object, with OTU, sample and taxa tables, after aggregation, for use with Shiny-Phyloseq (subfolder phyloseq)"
                     )
                   ),
                   tags$li(
                     tags$span(
                       "con: OTU and sample tables, after agregation, in a format suitable for use with the CoNet app (subfolder conet)"
                     )
                   ),
                   tags$li(
                     tags$span(
                       "gml: the bipartite sample and taxa network in .gml format, can be imported in Cytoscape and Gephi (subfolder gml)"
                     )
                   ),
                   tags$li(tags$span(
                     "ref: references for your selection (subfolder ref)"
                   ))
                 ),
                 style = "font-size:80%"
               ),
               p(
                 "A notification will briefly appear at the bottom right every time you export a file.",
                 "Please be patient: exporting aggregated tables may take quite a while for large selections.",
                 style = "font-size:80%"
               )
             ),
             # end sidebar panel Export tab

             # Main panel for displaying outputs ----
             mainPanel(
               # Output: Tabset w/ tables and plots ----
               tabsetPanel(
                 type = "tabs",
                 tabPanel(
                   "Summary",
                   br(),
                   h4("Summary info on your selection (before aggregation)"),
                   div(DT::dataTableOutput("mysummary"), style = "font-size:80%"),
                   hr(),
                   textOutput("agg_info")
                 ),
                 # end tab panel Summary
                 
                 tabPanel(
                   "OTU table",
                   h4("Your OTU table."),
                   p(
                     "Slider taxa2display lets you select how many taxa (after aggregation), ",
                     "to display. The top n taxa (n≤100, in decreasing order of relative abundance ",
                     "in your selection) are shown. No more than y taxa, were y is the number of ",
                     "taxa in your dataset after aggregation can be shown.",
                     style = "font-size:80%"
                   ),
                   p(
                     "Use the action button to save the table in the <table_plots> folder.",
                     style = "font-size:80%"
                   ),
                   hr(),
                   sliderInput(
                     "OTU",
                     label = h5("taxa2display"),
                     min = 0,
                     max = 100,
                     value = 20,
                     step = 1
                   ),
                   div(DT::dataTableOutput("OTU_table_agg"), style = "font-size:80%"),
                   hr(),
                   actionButton("save_OTU", "save"),
                   # diagnostic
                   # if (debug_mode) {verbatimTextOutput("sample_sums")},
                   if (debug_mode) {DT::dataTableOutput("OTU_test")}
                   
                 ),
                 # end tab panel OTU table  
                 
                 tabPanel(
                   "Prevalence table",
                   h4("Your prevalence and abundance table."),
                   p("The table is sorted by decreasing prevalence and abundance."),
                   div(DT::dataTableOutput("prev_ab_df"), style = "font-size:80%"),
                   hr(),
                   actionButton("save_prev_tab", "save")
                 ),
                 # end tab panel prevalence table
                 
                 tabPanel(
                   "Prev & ab plot",
                   fluidRow(
                     h4("Your prevalence and abundance graph."),
                     p(
                       "Slider <min.prev.> lets you select the fraction of samples in which a taxon ",
                       "needs to be present to be retained.  Slider <ab.tresh.> lets you select the ",
                       "abundance threshold: only OTU with maximum relative abundance larger than ",
                       "the treshold will be retained. Both filters affect the bar and the boxplot.",
                       style = "font-size:80%"
                     ),
                     p(
                       "Use the action button to save the graph in the <table_plots> folder. ",
                       "The graph will be saved in .jpg format, 150 dpi resolution.",
                       style = "font-size:80%"
                     ),
                     hr()
                   ),
                   # end first row tab prev ab plot
                   
                   fluidRow(
                     column(6,
                            sliderInput(
                              "prev",
                              label = h5("min.prev."),
                              min = 0,
                              max = 0.5,
                              value = minprev,
                              step = 0.01
                            )
                     ),
                     column(6,
                            sliderInput(
                              "minab",
                              label = h5("ab.tresh."),
                              min = 0,
                              max = 0.1,
                              value = minab,
                              step = 0.001
                            )
                     )
                   ),
                   # end second row tab prev ab plot
                   
                   fluidRow(
                     plotOutput("prev_ab_plot"),
                     hr(),
                     actionButton("save_prev_graph", "save"),
                     # diagnostic
                     if (debug_mode) {
                       verbatimTextOutput("pass_filter_text")
                     },
                     if (debug_mode) {
                       verbatimTextOutput("slider_prev_value")
                     },
                     if (debug_mode) {
                       verbatimTextOutput("slider_ab_value")
                     }
                   )
                   # end third row (graph and diagnostic) tab prev ab plot
                   
                 ),
                 # end tab prev ab plot
                 
                 tabPanel(
                   "Bar plot",
                   fluidRow(
                     h4("Your abundance plot."),
                     p(
                       "Menu <sample agg.> lets you select the sample aggregation level. ",
                       "Menu <taxa agg.> lets you select the taxa aggregation level. ",
                       "Both are affected by your choices in the <Aggregate> tab. ",
                       "Prevalence and abundance filters in the Prev & ab plot tab affect the barplot. ",
                       "No more than 25 taxa and 25 samples/samples groups are allowed. ",
                       style = "font-size:80%"
                     ),
                     p(
                       "Use the action button at the bottom to save the graph in the <table_plots> folder. ",
                       "The graph will be saved in .jpg format, 150 dpi resolution.",
                       style = "font-size:80%"
                     ),
                     hr()
                   ),
                   # end fist row (instructions) tab bar plot
                   
                   fluidRow(
                     column(6, uiOutput('sample_agg_lvl')),
                     column(6, uiOutput('taxa_agg_lvl'))
                   ),
                   # end second row (menus) tab bar plot
                   
                   fluidRow(
                     if (debug_mode) {
                       verbatimTextOutput("sample_menu_items")
                     },
                     if (debug_mode) {
                       verbatimTextOutput("taxa_menu_items")
                     }
                   ),
                   # end third row (diagnostic) tab bar plot
                   
                   fluidRow(
                     plotOutput("barplot"),
                     hr(),
                     actionButton("save_bar_plot", "save")
                   ),
                   # end fourth row tab bar plot                   
                   
                   fluidRow(
                     if(debug_mode) {
                       DT::dataTableOutput("barplot_table")
                     }
                   )
                   # end fifth row (diagnostic table) tab bar plot     
                 ),
                 # end tab bar plot
                 
                 tabPanel("Box plot",
                          fluidRow(
                            h4("Your (musical) taxa box plot."),
                            p(
                              "This plot is affected by choices in the filter, aggregate and bar plot tabs. ",
                              "Select box or violin from the <plot type> menu to produce either a box or a violin plot. ",
                              "Use the dropdown menus to choose the taxa aggregation level (first) and the taxa to show (last). ",
                              "Note that <all> taxa are shown by default; if you want a smaller selection of taxa, you must ",
                              "deselect <all> and manually select the taxa you want (at least 2).",
                              style = "font-size:80%"
                            ),
                            p(
                              "Use the action button at the bottom to save the graph in the <table_plots> folder. ",
                              "The graph will be saved in .jpg format, 150 dpi resolution.",
                              style = "font-size:80%"
                            ),
                            hr()
                          ),
                          # end first row (instructions) tab box plot
                          fluidRow(
                            column(4, uiOutput('taxa_agg_bxp')),
                            column(4, uiOutput('sel_taxa')),
                            column(4, 
                                   selectInput("g_type", "plot type", c("box", "violin"),
                                               selected = "box", multiple = F)
                                   )
                          ),
                          # end second tow, menus, tab box plot
                          fluidRow(
                          plotOutput("topxboxplot"),
                          hr(),
                          actionButton("save_box_plot", "save")
                          ),
                          # end second row (plot & save button) tab box plot
                          fluidRow(
                            if(debug_mode) {
                              DT::dataTableOutput("boxplot_table")
                            }
                          )
                          # end third row (diagnostic) tab box plot
                          
                 )
                 # end tab box plot
                 
               )
               # end tabset panel, main panel, Export tab
               
             )
             # end main panel, Export tab
             
           )
           # end sidebar layout, Export tab
  ),
  # end Export tab
  
  # About tab ---------------------------------------------------------------
  tabPanel("About",
           titlePanel("About"),
           sidebarLayout(
             sidebarPanel(
               h4("About FoodMicrobionet"),
               p(
                 "FoodMicrobionet is a database of studies on the
                                  composition of the bacterial microbiota of food obtained by
                                  amplicon targeted (16S RNA gene or 16S RNA) high throughput
                                  sequencing."
               ),
               p("Data are organized in tables:"),
               tags$div(tags$ul(
                 tags$li(tags$span("Studies: metadata on studies")),
                 tags$li(tags$span(
                   "Samples: metadata on samples for each study"
                 )),
                 tags$li(tags$span("Edges: taxa abundance for each sample")),
                 tags$li(tags$span("Taxa: metadata on taxa"))
               )),
               p(
                 "To learn more on FoodMicrobionet and its development visit ",
                 span(
                   a("our homepage. ",
                     href = "http://www.foodmicrobionet.org",
                     target = "_blank")
                 ),
                 "To cite FoodMicrobionet and see related articles ",
                 span(
                   a("go here.",
                     href = "https://doi.org/10.1016/j.ijfoodmicro.2015.12.001",
                     target = "_blank")
                 )
               ),
               h4("About this app"),
               p(FMBN$app_text),
               tags$div(tags$ul(
                 tags$li(tags$span("explore studies and samples;")),
                 tags$li(tags$span("extract/filter groups of samples;")),
                 tags$li(
                   tags$span("(optionally) perform aggregation of samples or taxa;")
                 ),
                 tags$li(tags$span("(locally) save data for further analysis;")),
                 tags$li(tags$span(
                   "generate and save some summary tables and graphs."
                 ))
               )),
               h4("Credits and copyright"),
               p(
                 "The function for installing/loading packages is taken verbatim from ",
                 span(
                   a(
                     "https://f1000research.com/articles/5-1492/v2",
                     href = "https://f1000research.com/articles/5-1492/v2",
                     target = "_blank"
                   )
                 ),
                 "."
               ),
               p(
                 "A few other pieces of code were adapted from ",
                 span(
                   a(
                     "https://tinyurl.com/tng8qgw",
                     href = "https://tinyurl.com/tng8qgw",
                     target = "_blank"
                   )
                 ),
                 "."
               ),
               h5(FMBN$copyright_text),
               p(
                 "Permission is hereby granted, free of charge, to any person obtaining a copy ",
                 "of this software and associated documentation files (the \"Software\"), to deal ",
                 "in the Software without restriction, including without limitation the rights ",
                 "to use, copy, modify, merge, publish, distribute, sublicense, and/or sell ",
                 "copies of the Software, and to permit persons to whom the Software is ",
                 "furnished to do so, subject to the following conditions:"
               ),
               tags$div(tags$ul(
                 tags$li(
                   tags$span(
                     "The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software."
                   )
                 ),
                 tags$li(
                   tags$span(
                     "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR ",
                     "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, ",
                     "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE ",
                     "AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER ",
                     "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, ",
                     "OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE ",
                     "SOFTWARE."
                   )
                 )
               )
               )
             ),
             # end sidebar panel, About tab
             
             mainPanel(
               p(
                 "A legacy image from FoodMicrobionet 3.1, generated with ",
                 span(a(
                   "Gephi.", href = "https://gephi.org",
                   target = "_blank"
                 )
                 )
               ),
               img(
                 src = "FMBN3_1_small.jpg",
                 height = 600
               )
             )
             # end main panel About tab
             
           ) 
           # end sidebar layout About tab
           
  ) 
  # end About tab
  
) 
  # end UI function
  
  
  # define server -----------------------------------------------------------
  
  server <- function(input, output, session) {
    # non reactive outputs -----------------------------------------------------
    # View tab, studies -------------------------------------------------------
    output$studies_table <- DT::renderDataTable(
      studies_disp,
      rownames = F,
      escape = c(8, 10),
      colnames = c(
        "study" = "studyId",
        "version" = "FMBN_version",
        "SRA/ENA" = "SRA_ENA",
        "reference" = "ref_short",
        "DOI link" = "DOIlink",
        "description" = "short_descr"
      ),
      options = list(
        pageLength = 3,
        lengthMenu = c(3, 5, 10, 15),
        columnDefs = list(list(
          targets = c(4, 9),
          render = JS(
            "function(data, type, row, meta) {",
            "return type === 'display' && data.length > 80 ?",
            "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
            "}"
          )
        ))
      )
    )
    
    # samples table, view tab -------------------------------------------------
    # Render selectInput for study, View tab
    
    output$choose_study <- renderUI({
      study_names <- as.list(studies$studyId)
      selectInput(
        "view_study",
        label = h5("Select one or more studies to explore samples"),
        choices = study_names,
        selected = "ST1",
        multiple = TRUE
      )
    })
    
    # Reactively subset so that only the selected study is viewed
    
    f_vsamples <- reactive({
      subset(
        samples,
        studyId %in% input$view_study,
        select = c(
          label_2,
          s_type,
          n_reads2,
          foodId,
          llabel,
          L1,
          L6,
          description
        )
      )
    })
    
    output$vsamples_table <- DT::renderDataTable(
      f_vsamples(),
      rownames = F,
      escape = T,
      colnames = c(
        "label" = "label_2",
        "sample type" = "s_type",
        "reads" = "n_reads2",
        "code" = "foodId",
        "ext. code" = "llabel",
        "food group" = "L1",
        "food subgroup" = "L6",
        "description" = "description"
      ),
      options = list(
        pageLength = 5,
        lengthMenu = c(5, 10, 20, 25),
        columnDefs = list(# limits display length of column 5 = description
          list(
            targets = 6,
            render = JS(
              "function(data, type, row, meta) {",
              "return type === 'display' && data.length > 30 ?",
              "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
              "}"
            )
          ))
        
      )
    )
    
    # functions for the filter tab --------------------------------------------
    
    # gets the list of studies to be selected
    study_sel <- reactive({
      if (is.null(input$select_study) | input$study_ch_box == F) {
        studies$studyId
      } else {
        input$select_study
      }
    },
    label = "study_sel_reactive")
    
    # reactive for the dynamic food group list
    food_group_names <- reactive({
      nfstudy <- subset(samples, studyId %in% study_sel())
      unique(nfstudy$L1)
    })
    
    # renders the food group drop down
    output$choose_food_group <- renderUI({
      selectInput(
        "select_L1",
        label = "Select one or more food groups",
        choices = as.list(food_group_names()),
        multiple = TRUE
      )
    })
    
    # creates the reactive which responds to select_L1
    food_group_sel <- reactive({
      if (is.null(input$select_L1) | input$food_group_ch_box == F) {
        # using dplyr::pull to pull a vector
        food_group_names()
      } else {
        input$select_L1
      }
    },
    label = "food_group_sel_reactive")
    
    # reactive, for the dynamic foodId list
    
    foodId_names <- reactive({
      nfstudy <- subset(samples, studyId %in% study_sel())
      nffgroup <- subset(nfstudy, L1 %in% food_group_sel())
      unique(nffgroup$foodId)
    })
    
    # renders the foodId drop down
    output$choose_foodId <- renderUI({
      selectInput(
        "select_foodId",
        label = "Select one or more food codes",
        choices = as.list(foodId_names()),
        multiple = TRUE
      )
    })
    
    # creates the reactive which responds to select_foodId
    foodId_sel <- reactive({
      if (is.null(input$select_foodId) | input$foodId_ch_box == F) {
        foodId_names()
      } else {
        input$select_foodId
      }
    },
    label = "foodId_sel_reactive")
    
    # creates the reactive to display the table
    f_fsamples <- reactive({
      fstudy <- subset(ssamples, studyId %in% study_sel())
      ffgroup <- subset(fstudy, L1 %in% food_group_sel())
      ffoodId <- subset(ffgroup, foodId %in% foodId_sel())
      ffinal <- ffoodId
      ffinal <- mutate(ffinal, SRA_run = ifelse(
        is.na(SRA_run),
        "not available",
        paste0(
          "<a href='",
          "https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=",
          SRA_run,
          "' target='_blank'>",
          SRA_run,
          "</a>"
        )
      ))
      return(ffinal)
    })
    
    # renders the table for the filter tab
    output$ssamplestable <- DT::renderDataTable(
      f_fsamples(),
      filter = "bottom",
      rownames = F,
      escape = 21,
      colnames = c(
        "study" = "studyId",
        "label" = "label_2",
        "sample type" = "s_type",
        "reads" = "n_reads2",
        "issues" = "n_issues",
        "code" = "foodId",
        "ext. code" = "llabel",
        "L1 food group" = "L1",
        "L6 food subgroup" = "L6",
        "description" = "description",
        "DNA/cDNA" = "target1",
        "region" = "target2",
        "run" = "SRA_run"
      ),
      options = list(
        pageLength = 10,
        lengthMenu = c(5, 10, 20, 30, 40),
        scrollX = T,
        columnDefs = list(
          # limits the variables which are visible
          list(
            visible = FALSE,
            targets = c(1, 2, 4, 5, 7, 13, 20, 21)
          ),
          # defines which columns are searchable
          list(
            targets = c(0, 2, 
                        3, 10, 11, 12, 14, 22),
            searchable = FALSE
          ),
          # limits number of characters shown, more on hover
          list(
            targets = c(11, 12, 13, 15),
            render = JS(
              "function(data, type, row, meta) {",
              "return type === 'display' && data.length > 20 ?",
              "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
              "}"
            )
          )
        )
        
      )
    )
    
    # diagnostic, activate by setting debug_mode <- T at line 92
    if (debug_mode) {
      output$filter_ssampletable_rows <-
        renderPrint(input$ssamplestable_rows_all)
      samples_label1 <- eventReactive(input$ssamplestable_rows_all, {
        # extract samples based on row numbers of the filtered table
        samples_filt1 <-
          pull(f_fsamples()[input$ssamplestable_rows_all, 3])
        return(samples_filt1)
      })
      output$samples_label <- renderPrint(samples_label1())
      samples_id_a <- eventReactive(input$ssamplestable_rows_all, {
        # extract samples based on row numbers of the filtered table
        samples_filt2 <-
          pull(f_fsamples()[input$ssamplestable_rows_all, 2])
        return(samples_filt2)
      })
      output$samples_ids <- renderPrint(samples_id_a())
    }
    
    # functions for the aggregate tab -----------------------------------------
    # diagnostic, activate by setting debug_mode <- T in options
    if (debug_mode) {
      # row numbers from the filtered table (in the filter tab)
      row_numbers <- eventReactive(input$ssamplestable_rows_all, {
        pull(f_fsamples()[input$ssamplestable_rows_all, 2])
      })
      output$sel_row_numbers <- renderPrint(row_numbers())
      sel_samples_slice <- eventReactive(input$ssamplestable_rows_all, {
        rows_to_slice <- pull(f_fsamples()[input$ssamplestable_rows_all, 2])
        # a subset from samples
        return(dplyr::slice(dplyr::select(samples, studyId:label_1),
                            rows_to_slice))
      })
      output$agg_samples_mt <-
        DT::renderDataTable(sel_samples_slice())
      # end diagnostic
    }
    
    FMBNfilt <- eventReactive(
      c(input$ssamplestable_rows_all,
        input$study_ch_box,
        input$food_group_ch_box,
        input$foodId_ch_box
        ), 
      {
      if(any(c(input$study_ch_box, input$food_group_ch_box, input$foodId_ch_box))){
      # extract samples based on row numbers of the filtered table
      # this only works if all sampleIds are consecutive and match row numbers
      # a slower alternative would be selection based on sample labels (label_1)
      # using perhaps the %in% operator, but performs much worse in benchmarking
      samples_filt <- samples %>%
        dplyr::slice(pull(f_fsamples()[input$ssamplestable_rows_all, 2]))
      # get studies
      studies_filt <- studies %>%
        dplyr::filter(studyId %in% pull(distinct(
          dplyr::select(samples_filt, studyId), studyId
        )))
      # get references
      reference_filt <- dplyr::bind_rows(references,
                                         dplyr::select(studies_filt,
                                                       studyId, ref_complete, DOI))
      # get edges
      edges_filt <-
        edges %>% dplyr::filter(sampleId %in% samples_filt$sampleId)
      
      # get relative abundance of taxa
      taxon_ab <- edges_filt %>% group_by(taxonId) %>%
        summarise(ab_sum = sum(weight)) %>%
        ungroup() %>%
        mutate(ab_sum = ab_sum / sum(ab_sum))
      
      # get taxa from edges
      taxonId_filt <-
        edges_filt %>% dplyr::distinct(taxonId) %>% arrange(.) %>% pull(.)
      
      # filter taxa
      taxa_filt <-
        taxa %>% dplyr::filter(taxonId %in% taxonId_filt) %>%
        left_join(., taxon_ab)
      
      # building a list with essential elements ---------------------------------
      # initially sample aggregation is sample and taxonomic aggregation is species
      
      FMBN_list <- list(
        studies = studies_filt,
        samples = samples_filt,
        edges = edges_filt,
        taxa = taxa_filt,
        references = reference_filt,
        version = paste0(
          "This dataset includes ",
          nrow(samples_filt),
          " samples from FoodMicrobionet, extracted on ",
          today(),
          ". ",
          version
        ),
        sample_Ids = pull(f_fsamples()[input$ssamplestable_rows_all, 2]),
        sample_agg = "sample",
        tax_agg = "species"
      )
      
      } else {
        FMBN_list <- NULL
        showNotification(
          "ERROR you must select at least one of the checkboxes in the filter tab",
          type = "message",
          duration = 10,
          closeButton = TRUE
        )
      }
        
        return(FMBN_list)
    })
    
    # building the reactive aggregate list -------------------------------------
    # will be used in any case to export data
    
    FMBNexp <- eventReactive(
      c(input$ssamplestable_rows_all,
        input$tax_agg_level,
        input$sample_agg_level,
        input$study_ch_box,
        input$food_group_ch_box,
        input$foodId_ch_box
      ),
      {
        if(any(c(input$study_ch_box, input$food_group_ch_box, input$foodId_ch_box))) {
          # conditions for sample aggregation
          col_agg_var <- switch(input$sample_agg_level,
                                "sample" = "label_2",
                                "exp. code" = "llabel")
          # handle the case when there is only one exp. code in the selection
          no_agg_flag <- F
          actual_sample_agg <- input$sample_agg_level
          do_no_agg <- (col_agg_var == "llabel" & n_distinct(FMBNfilt()[[2]]$llabel)<2)
          
          if(do_no_agg){
            col_agg_var <- "label_2"
            no_agg_f <- T
            actual_sample_agg <- "sample"
            showNotification(
              "too few exp. code categories, no aggregation possible",
              type = "message",
              duration = 10,
              closeButton = TRUE
            )
          }
          
          taxa_redux <- switch(
            input$tax_agg_level,
            "species" = {
              FMBNfilt()[[4]] %>%
                mutate(a_label = label, t_label = genus)
            },
            "genus" = {
              FMBNfilt()[[4]] %>%
                dplyr::select(taxonId:taxonomy, idelevel) %>%
                mutate(species = NA)
            },
            "family" = {
              FMBNfilt()[[4]] %>%
                dplyr::select(taxonId:taxonomy, idelevel) %>%
                mutate(genus = NA, species = NA)
            },
            "class" = {
              FMBNfilt()[[4]] %>%
                dplyr::select(taxonId:taxonomy, idelevel) %>%
                mutate(
                  order = NA,
                  family = NA,
                  genus = NA,
                  species = NA
                )
            }
          )
          if (input$tax_agg_level != "species") {
            taxa_redux <- switch(
              input$tax_agg_level,
              "genus" = {
                taxa_redux %>%
                  mutate(
                    id_L6 = str_sub(id_L6, 1, str_locate(id_L6, ";s_")[, 1] + 3),
                    taxonomy = str_sub(taxonomy, 1, str_locate(taxonomy, "; s_")[, 1] +
                                         4),
                    t_label = genus,
                    a_label = ifelse(is.na(t_label), label, t_label)
                  )
              },
              "family" = {
                taxa_redux %>%
                  mutate(
                    id_L6 = paste0(str_sub(
                      id_L6, 1, str_locate(id_L6, ";g_")[, 1] + 3
                    ), ";s__"),
                    taxonomy = paste0(str_sub(
                      taxonomy, 1, str_locate(taxonomy, "; g_")[, 1] + 4
                    ), "; s__"),
                    t_label = family,
                    a_label = ifelse(is.na(t_label), label, t_label)
                  )
              },
              "class" = {
                taxa_redux %>%
                  mutate(
                    id_L6 = paste0(str_sub(
                      id_L6, 1, str_locate(id_L6, ";f_")[, 1] + 3
                    ), ";g__;s__"),
                    taxonomy = paste0(str_sub(
                      taxonomy, 1, str_locate(taxonomy, "; f_")[, 1] + 4
                    ), "; g__; s__"),
                    t_label = class,
                    a_label = ifelse(is.na(t_label), label, t_label)
                  )
              }
            )
          }
          # adding info to the edges
          edges_exp <- left_join(
            FMBNfilt()[[3]],
            dplyr::select(taxa_redux,-id_L6,-idelevel,-label,-taxonomy,-t_label)
          )
          edges_exp <- left_join(edges_exp,
                                 dplyr::select(FMBNfilt()[[2]], sampleId, label_2, llabel, foodId))
          
          # cast
          
          OTUtableagg <-
            dcast(
              edges_exp,
              as.formula(paste("a_label", col_agg_var, sep = "~")),
              sum,
              value.var = "weight",
              drop = FALSE,
              fill = 0
            )
          
          # transforming in proportions
          sumcheck <- colSums(OTUtableagg[2:ncol(OTUtableagg)])
          tOTUm <- t(OTUtableagg[, 2:ncol(OTUtableagg)]) / sumcheck
          colnames(tOTUm) <- OTUtableagg$a_label
          OTUtableagg[, 2:ncol(OTUtableagg)] <- t(tOTUm)
          taxa_list <-
            taxa_redux %>% distinct(a_label, .keep_all = T) %>%
            dplyr::select(-taxonId,-label,-t_label) %>% 
            mutate(idelevel = case_when(
              is.na(domain) | is.na(phylum) ~ "domain",
              is.na(class) ~ "phylum",
              is.na(order) ~ "class",
              is.na(family) ~ "order",
              is.na(genus)  ~ "family",
              is.na(species) ~ "genus",
              TRUE ~ "species"
            ))
          
          # getting nreads
          nreads <- switch(
            col_agg_var,
            "label_2" = FMBNfilt()[[2]] %>%
              dplyr::select(label = label_2, n_reads2) %>% arrange(label),
            "llabel" = FMBNfilt()[[2]] %>%
              dplyr::select(label = llabel, n_reads2) %>% arrange(label) %>%
              mutate(n_reads2 = median(n_reads2)) %>%
              group_by(label) %>%
              summarise(n_reads2 = median(n_reads2) * n())
          )
          # transforming tOTUm
          tOTUm <- round(tOTUm * nreads$n_reads2)
          # prepping the sample table
          sample_table <- switch(col_agg_var,
                                 "label_2" = {
                                   FMBNfilt()[[2]] %>%
                                     dplyr::select(studyId, label = label_2, llabel, s_type,
                                                   n_reads2:SRA_run)
                                 },
                                 "llabel" = {
                                   FMBNfilt()[[2]] %>%
                                     dplyr::select(-n_reads2) %>%
                                     mutate(
                                       studyId = NA,
                                       label = llabel,
                                       description = " ",
                                       target1 = NA,
                                       target2 = NA
                                     ) %>%
                                     distinct(label, .keep_all = TRUE) %>%
                                     left_join(., nreads, by = "label") %>%
                                     dplyr::select(studyId, label, llabel, s_type, n_reads2,
                                                   foodId:SRA_run)
                                 })
          
          sample_table <- as.data.frame(sample_table)
          rownames(sample_table) <- sample_table$label
          
          # making an object for import in Gephi  or Cytoscape ---------------------------
          
          edges <-
            OTUtableagg %>% tidyr::gather(key = "Source", value = "weight",-a_label) %>%
            dplyr::mutate(weight = weight * 100) %>%
            dplyr::select(Source, Target = a_label, weight) %>%
            dplyr::filter(weight > 0)
          
          sample_nodes <- switch(col_agg_var,
                                 "label_2" = {
                                   sample_table %>%
                                     dplyr::select(
                                       label,
                                       node_type = s_type,
                                       L1:L6,
                                       info1 = llabel,
                                       info2 = foodId,
                                       info3 = nature,
                                       info4 = process,
                                       info5 = spoilage,
                                       info6 = description,
                                       info7 = target1,
                                       info8 = target2
                                     )
                                 },
                                 "llabel" = {
                                   sample_table %>%
                                     mutate(
                                       info1 = label,
                                       info6 = as.character(NA),
                                       info7 = as.character(NA),
                                       info8 = as.character(NA)
                                     ) %>%
                                     dplyr::select(
                                       label,
                                       node_type = s_type,
                                       L1:L6,
                                       info1,
                                       info2 = foodId,
                                       info3 = nature,
                                       info4 = process,
                                       info5 = spoilage,
                                       info6:info8
                                     )
                                   
                                 })
          
          OTU_nodes <- taxa_list %>%
            mutate(
              node_type = "OTU",
              L1 = phylum,
              L4 = class,
              L6 = family,
              info1 = genus,
              info2 = species,
              info3 = NA,
              info4 = NA,
              info5 = NA,
              info6 = id_L6,
              info7 = NA,
              info8 = NA
            ) %>%
            dplyr::select(label = a_label, node_type, L1:info8)
          
          nodes <- bind_rows(sample_nodes, OTU_nodes)
          
          # make a igraph object
          FMBNigraph <-
            graph_from_data_frame(edges, directed = F, vertices = nodes)
          
          # adding abundance to taxa_list
          taxa_ab <- edges %>% group_by(Target) %>%
            summarise(ab_sum = sum(weight)) %>%
            ungroup() %>%
            mutate(ab_sum = ab_sum / sum(ab_sum)) %>%
            dplyr::rename(a_label = Target)
          taxa_list <- left_join(taxa_list, taxa_ab)
          
          FMBN_agg <- list(
            studies = FMBNfilt()[[1]],
            OTU_table = tOTUm,
            OTU_table_relf = OTUtableagg,
            sample_metadata = sample_table,
            taxa_metadata = taxa_list %>% dplyr::select(label = a_label,
                                                        domain:taxonomy, ab_sum, idelevel),
            edge_table = edges,
            node_table = nodes,
            i_graph = FMBNigraph,
            references = FMBNfilt()[[5]],
            version = FMBNfilt()[[6]],
            sample_agg = actual_sample_agg,
            tax_agg = input$tax_agg_level,
            diagnostic_f = do_no_agg
          )
        } else {
          FMBN_agg <- NULL
          showNotification(
            "ERROR you must select at least one of the checkboxes in the filter tab",
            type = "message",
            duration = 10,
            closeButton = TRUE
          )                         
        }
        return(FMBN_agg)
      })

    if(debug_mode) output$no_agg <- renderPrint(FMBNexp()[[13]])
    
    # the outputs, aggregate tab ----------------------------------------------
    
    # data on selection, aggregate tab
    output$subset_info <- renderText(FMBNfilt()[[6]])
    
    # study table, aggregate tab, not affected by aggregation
    output$f_studies_table <- DT::renderDataTable({
      FMBNfilt()[[1]] %>%
        mutate(
          target = paste(target, region),
          pipeline = str_c(
            bioinf_software,
            OTU_picking,
            assign_tax_method,
            tax_database,
            sep = "; "
          ),
          SRA_ENA = ifelse(
            is.na(Seq_accn),
            "not available",
            paste0(
              "<a href='",
              Seq_accn_link,
              "' target='_blank'>",
              Seq_accn,
              "</a>"
            )
          ),
          DOIlink = ifelse(
            DOI == "unpublished data",
            "unpublished data",
            paste0("<a href='",
                   DOI_link,
                   "' target='_blank'>",
                   DOI, "</a>")
          )
        ) %>%
        dplyr::select(
          studyId,
          FMBN_version,
          target,
          platform,
          pipeline,
          samples,
          SRA_ENA,
          ref_short,
          DOIlink,
          short_descr
        )
    },
    rownames = F,
    escape = c(8, 10),
    colnames = c(
      "study" = "studyId",
      "version" = "FMBN_version",
      "SRA/ENA" = "SRA_ENA",
      "reference" = "ref_short",
      "Seq. accn. link" = "SRA_ENA",
      "DOI link" = "DOIlink",
      "description" = "short_descr"
    ),
    options = list(
      pageLength = 3,
      lengthMenu = c(3, 5, 10),
      columnDefs = list(list(
        targets = c(4, 9),
        render = JS(
          "function(data, type, row, meta) {",
          "return type === 'display' && data.length > 50 ?",
          "'<span title=\"' + data + '\">' + data.substr(0, 50) + '...</span>' : data;",
          "}"
        )
      ))
    ))
    # sample table, aggregate tab, affected by sample_agg_level and show_what
    
    output$f_samples_table <- DT::renderDataTable(
      if (input$sample_agg_level == "samples" |
          input$show_what == "original") {
        dplyr::select(
          FMBNfilt()[[2]],
          studyId,
          label = label_2,
          s_type,
          n_reads2,
          foodId,
          llabel,
          L1,
          L6,
          description,
          nature:target2
        )
      } else {
        dplyr::select(
          FMBNexp()[[4]],
          studyId,
          label,
          s_type,
          n_reads2,
          foodId,
          llabel,
          L1,
          L6,
          description,
          nature:target2
        )
      },
      rownames = F,
      escape = T,
      colnames = c(
        "study" = "studyId",
        "label" = "label",
        "sample type" = "s_type",
        "reads" = "n_reads2",
        "code" = "foodId",
        "exp. code" = "llabel",
        "L1 food group" = "L1",
        "L6 food subgroup" = "L6",
        "description" = "description",
        "DNA/cDNA" = "target1",
        "region" = "target2"
      ),
      options = list(
        pageLength = 3,
        lengthMenu = c(3, 5, 10, 20, 25),
        columnDefs = list(# limits display length of columns 7-9
          list(
            targets = c(7, 8, 9),
            render = JS(
              "function(data, type, row, meta) {",
              "return type === 'display' && data.length > 20 ?",
              "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
              "}"
            )
          ))
        
      )
    )
    
    # taxa table, aggregate tab, affected by sample_agg_level and show_what
    
    output$f_taxa_table <- DT::renderDataTable(
      if (input$tax_agg_level == "species" |
          input$show_what == "original") {
        FMBNfilt()[[4]] %>%
          dplyr::arrange(-ab_sum) %>%
          dplyr::select(label,
                        ab_sum,
                        domain:species,
                        NCBI_outlink,
                        BacterioNet_outlink,
                        Flori_habitat,
                        Flori_pheno,
                        Flori_use,
                        idelevel) %>%
          dplyr::mutate(
            ab_sum = round(ab_sum, digits = 4),
            NCBI_outlink = str_c(
              "<a href='",
              NCBI_outlink,
              "' target='_blank'>",
              "go!",
              "</a>",
              sep=""
            ),
            BacterioNet_outlink = str_c(
              "<a href='",
              BacterioNet_outlink,
              "' target='_blank'>",
              "go!",
              "</a>",
              sep=""
            ),
            Flori_habitat = if_else(
              idelevel %in% c("domain", "phylum", "class", "order"),
              "no search path",
              str_c(
                "<a href='",
                Flori_habitat,
                "' target='_blank'>",
                "go!",
                "</a>",
                sep =""
              )
            ),
            Flori_pheno = if_else(
              idelevel %in% c("domain", "phylum", "class", "order"),
              "no search path",
              str_c(
                "<a href='",
                Flori_pheno,
                "' target='_blank'>",
                "go!",
                "</a>",
                sep=""
              )
            ),
            Flori_use = if_else(
              idelevel %in% c("domain", "phylum", "class", "order"),
              "no search path",
              str_c(
                "<a href='",
                Flori_use,
                "' target='_blank'>",
                "go!",
                "</a>",
                sep=""
              )
            )
          ) %>%
          select(-idelevel)
      } else {
        FMBNexp()[[5]] %>%
          dplyr::arrange(-ab_sum) %>%
          dplyr::select(label, ab_sum, domain:species, idelevel) %>%
          dplyr::mutate(
            ab_sum = round(ab_sum, digits = 4),
            NCBI_outlink = case_when(
              label == "Other" | label == "" ~ str_c("<a href='",
                                                     "http://www.foodmicrobionet.org/?p=235",
                                                     "' target='_blank'>",
                                                     "go!",
                                                     "</a>",
                                                     sep=""),
              str_detect(label, " \\(class\\)") ~ str_c("<a href='",
                                                        "http://www.ncbi.nlm.nih.gov/taxonomy/?term=", 
                                                        str_remove(label, pattern =" \\(class\\)"),
                                                        "' target='_blank'>",
                                                        "go!",
                                                        "</a>",
                                                        sep=""),
              TRUE ~ str_c("<a href='",
                           "http://www.ncbi.nlm.nih.gov/taxonomy/?term=", 
                           label, 
                           "' target='_blank'>",
                           "go!",
                           "</a>",
                           sep ="")
            ),
            BacterioNet_outlink = case_when(
              label == "Other" | label == "" ~ str_c("<a href='",
                                                    "http://www.foodmicrobionet.org/?p=235",
                                                    "' target='_blank'>",
                                                    "go!",
                                                    "</a>",
                                                    sep=""),
              str_detect(label, " \\(class\\)") ~ str_c("<a href='",
                                                        "https://lpsn.dsmz.de/", 
                                                        idelevel, 
                                                        "/", 
                                                        str_remove(string = label, pattern =" \\(class\\)"),
                                                        "' target='_blank'>",
                                                        "go!",
                                                        "</a>",
                                                        sep=""),
              TRUE ~ str_c("<a href='",
                           "https://lpsn.dsmz.de/", 
                           idelevel, 
                           "/", 
                           label, 
                           "' target='_blank'>",
                           "go!",
                           "</a>",
                           sep ="")
            ),
            Flori_habitat = case_when(
              label == "Other" | label == "" ~ str_c("<a href='",
                                                     "http://www.foodmicrobionet.org/?p=235",
                                                     "' target='_blank'>",
                                                     "go!",
                                                     "</a>",
                                                     sep=""),
              idelevel %in% c("domain", "phylum", "class", "order") ~ "no search path",
              TRUE ~ str_c("<a href='",
                           "http://migale.jouy.inra.fr/Florilege/#&searchByTaxon=", 
                           label, 
                           "' target='_blank'>",
                           "go!",
                           "</a>",
                           sep ="")
            ),
            Flori_pheno = case_when(
              label == "Other" | label == "" ~ str_c("<a href='",
                                                     "http://www.foodmicrobionet.org/?p=235",
                                                     "' target='_blank'>",
                                                     "go!",
                                                     "</a>",
                                                     sep=""),
              idelevel %in% c("domain", "phylum", "class", "order") ~ "no search path",
              TRUE ~ str_c("<a href='",
                           "http://migale.jouy.inra.fr/Florilege/#&searchByTaxonForPhenotype=", 
                           label, 
                           "' target='_blank'>",
                           "go!",
                           "</a>",
                           sep ="")
            ),
            Flori_use =case_when(
              label == "Other" | label == "" ~ str_c("<a href='",
                                                     "http://www.foodmicrobionet.org/?p=235",
                                                     "' target='_blank'>",
                                                     "go!",
                                                     "</a>",
                                                     sep=""),
              idelevel %in% c("domain", "phylum", "class", "order") ~ "no search path",
              TRUE ~ str_c("<a href='",
                           "http://migale.jouy.inra.fr/Florilege/#&searchByTaxonForUse=", 
                           label, 
                           "' target='_blank'>",
                           "go!",
                           "</a>",
                           sep ="")
            )
          ) %>%
          select(-idelevel)
      },
      rownames = F,
      escape = c(5, 6, 7, 8, 9),
      colnames = c(
        "label" = "label",
        "rel. abund." = "ab_sum",
        "domain" = "domain",
        "phylum" = "phylum",
        "class" = "class",
        "family" = "family",
        "genus" = "genus",
        "species" = "species",
        "NCBI" = "NCBI_outlink",
        "LPSN" = "BacterioNet_outlink",
        "Fl.habitat" = "Flori_habitat",
        "Fl.pheno" = "Flori_pheno",
        "Fl.use" = "Flori_use"
      ),
      options = list(
        pageLength = 3,
        lengthMenu = c(3, 5, 10, 20, 25)
      )
    )
    # reference list, aggregate tab
    output$f_refs <- DT::renderDataTable({
      FMBNfilt()[[5]] %>% mutate(DOIlink = ifelse(
        DOI == "unpublished data",
        "unpublished data",
        paste0(
          "<a href='",
          paste0("https://doi.org/", DOI),
          "' target='_blank'>",
          DOI,
          "</a>"
        )
      )) %>%
        dplyr::select(-DOI)
    },
    rownames = F,
    escape = c(2),
    colnames = c(
      "Set" = "studyId",
      "Reference." = "ref_complete",
      "DOI" = "DOIlink"
    ),
    options = list(
      pageLength = 5,
      lengthMenu = c(3, 5, 10, 20, 25)
    ))
    
    # functions, export tab ---------------------------------------------------
    
    # reactives for saving object from the sidebar pane -----------------------
    
    # save selection no aggregation
    observeEvent(input$miniFMBN,
                 {
                   saveRDS(FMBNfilt(),
                           file = file.path(
                             "output",
                             "minifmbn",
                             paste(input$fn_prefix, input$miniFMBN, ".RDS", sep =
                                     "")
                           ))
                   showNotification(
                     "urFMBN file saved to /output/minifmbn folder",
                     type = "message",
                     duration = 5,
                     closeButton = TRUE
                   )
                 })
    
    # save selection after aggregation
    observeEvent(input$aggFMBN,
                 {
                   saveRDS(FMBNexp(),
                           file = file.path(
                             "output",
                             "aggdata",
                             paste(input$fn_prefix, input$aggFMBN, "_agg.RDS", sep =
                                     "")
                           ))
                   showNotification(
                     "aggFMBN file saved to /output/aggdata folder",
                     type = "message",
                     duration = 5,
                     closeButton = TRUE
                   )
                 })
    
    # create and save save a phyloseq object after aggregation
    observeEvent(input$phyloseq,
                 {
                   taxtable <- FMBNexp()[[5]] %>% dplyr::select(label, domain:species)
                   taxtable <- as.data.frame(taxtable)
                   rownames(taxtable) <- taxtable$label
                   taxtable <- as.matrix(taxtable[, 2:ncol(taxtable)])
                   # create and a phyloseq class object
                   physeqdata <-
                     phyloseq(otu_table(t(FMBNexp()[[2]]), taxa_are_rows = T),
                              tax_table(taxtable),
                              sample_data(FMBNexp()[[4]]))
                   # phyloseq object
                   save(physeqdata,
                        file = file.path(
                          "output",
                          "phyloseq",
                          paste(input$fn_prefix,
                                input$phyloseq,
                                "_physeq.Rdata", sep =
                                  "")
                        ))
                   showNotification(
                     "phyloseq file saved to /output/phyloseq folder",
                     type = "message",
                     duration = 5,
                     closeButton = TRUE
                   )
                 })
    
    # create and save a OTU and sample table for CoNet
    observeEvent(input$conet,
                 {
                   OTU_table_CoNet <- as.data.frame(t(FMBNexp()[[2]])) %>%
                     tibble::rownames_to_column(var = "OTUID") %>%
                     inner_join(.,
                                dplyr::select(FMBNexp()[[5]], taxonomy, label),
                                by = c("OTUID" = "label"))
                   write_tsv(OTU_table_CoNet,
                             file.path(
                               "output",
                               "conet",
                               paste(input$fn_prefix,
                                     input$conet,
                                     "_conet_OTU.txt", sep =
                                       "")
                             ))
                   write_tsv(
                     remove_rownames(FMBNexp()[[4]]) %>%
                       dplyr::rename(SampleID = label) %>%
                       dplyr::select(SampleID, setdiff(names(FMBNexp(
                       )[[4]]), "label")),
                     file.path(
                       "output",
                       "conet",
                       paste(input$fn_prefix,
                             input$conet,
                             "_conet_samples.txt", sep =
                               "")
                     )
                   )
                   showNotification(
                     "OTU and sample tables for CoNet saved to /output/conet folder",
                     type = "message",
                     duration = 5,
                     closeButton = TRUE
                   )
                 })
    
    # create and save a .gml graph for use in Gephi and Cytoscape
    observeEvent(input$gml,
                 {
                   write.graph(FMBNexp()[[8]],
                               file.path(
                                 "output",
                                 "gml",
                                 paste(input$fn_prefix,
                                       input$gml,
                                       "_igraph.gml", sep =
                                         "")
                               ),
                               format = "gml")
                   showNotification(
                     ".gml network saved in /output/gml folder",
                     type = "message",
                     duration = 5,
                     closeButton = TRUE
                   )
                 })
    # export a tab delimited files with the references for your selection
    observeEvent(input$ref,
                 {
                   write_tsv(FMBNfilt()[[5]],
                             file.path(
                               "output",
                               "ref",
                               paste(input$fn_prefix, input$ref,
                                     "_references.txt", sep = "")
                             ))
                   showNotification(
                     "references saved in /output/ref folder",
                     type = "message",
                     duration = 5,
                     closeButton = TRUE
                   )
                 })
    
    # for debugging
    if (debug_mode)
      output$show_output <- renderText({
        input$show_summary
      })
    
    
    # building objects for the output tabs ------------------------------------
    # resetting some inputs in this pane if the filtered table change
    # NOTE rarefaction should be applied at this stage and a control should be added in the sidebar
    
    # the summary info (this observes several inputs)
    my_summaries <- eventReactive(
      c(
        input$show_summary,
        input$ssamplestable_rows_all,
        input$sample_agg_level,
        input$tax_agg_level
      ),
      {
        if (!input$show_summary)
        {
          NULL
        } else {
          # make the OTU table
          # remove Eukaryota, Chloroplast, Mitochondria, Other
          euk_labels <- FMBNexp()[[5]] %>%
            dplyr::filter(domain == "Eukaryota" | family == "Mitochondria" | class == "Chloroplast") %>%
            pull(label)
          remove_what <-
            c(euk_labels, "Other", "Chloroplast", "Mitochondria")
          to_remove <-
            which(colnames(FMBNexp()[[2]]) %in% remove_what)
          if(length(to_remove)>0){
            my_OTU_table <- FMBNexp()[[2]][,-to_remove]
          } else {
            my_OTU_table <- FMBNexp()[[2]]
          }
          my_OTU_table_relf <- my_OTU_table / rowSums(my_OTU_table)
          # make the prevalence table
          prevdf <- apply(
            X = my_OTU_table,
            MARGIN = 2,
            FUN = function(x) {
              sum(x > 0)
            }
          )
          min_rel_ab <- apply(
            X = my_OTU_table_relf,
            MARGIN = 2,
            FUN = function(x) {
              min(x)
            }
          )
          max_rel_ab <- apply(
            X = my_OTU_table_relf,
            MARGIN = 2,
            FUN = function(x) {
              max(x)
            }
          )
          rel_prev <- prevdf / ncol(FMBNexp()[[3]])
          # Add taxonomy and total read counts to this data.frame
          prevdf <- data.frame(
            Prev = prevdf,
            Rel_Prev = round(rel_prev, 3),
            Tot_Ab = colSums(my_OTU_table),
            Min_Rel_Ab = round(min_rel_ab, 4),
            Max_Rel_Ab = round(max_rel_ab, 4)
          )
          prevdf <- prevdf %>%
            rownames_to_column(var = "label") %>%
            mutate(Rel_Ab = round(Tot_Ab / sum(Tot_Ab), 4)) %>%
            left_join(., dplyr::select(FMBNexp()[[5]], label, phylum:species)) %>%
            dplyr::arrange(desc(Prev), desc(Rel_Ab))
          
          # the list returned by this eventReactive
          # 1. the summary table
          # 2. the summary text
          # 3. the OTU table, abs. abundance (number of sequences)
          # 4. taxa (with taxonomic information), ordered in decreasing order of abundance
          # 5. a prevalence and abundance table
          summary_info <- list(
            summary_table <-  data.frame(
              Summary = c(
                "Studies",
                "Samples",
                "Food groups",
                "Food codes",
                "exp. food codes",
                "taxa"
              ),
              number = c(
                nrow(FMBNfilt()[[1]]),
                nrow(FMBNfilt()[[2]]),
                n_distinct(FMBNfilt()[[2]]$L1),
                n_distinct(FMBNfilt()[[2]]$foodId),
                n_distinct(FMBNfilt()[[2]]$llabel),
                nrow(FMBNfilt()[[4]])
              )
            ),
            summary_text <- paste(
              "The sample aggregation level is ",
              FMBNexp()[[11]],
              ". ",
              "The aggregation level for taxa is ",
              input$tax_agg_level,
              ". There are ",
              nrow(FMBNexp()[[3]]),
              " taxa and ",
              ncol(FMBNexp()[[3]]),
              " samples/sample groups in your",
              " dataset after aggregation.",
              sep = ""
            ),
            OTU_table = my_OTU_table,
            f_top_taxa = FMBNexp()[[5]] %>%
              dplyr::filter(!(label %in% remove_what)) %>%
              dplyr::arrange(desc(ab_sum)),
            prev_table <- prevdf
          )
          return(summary_info)
        }
      }
    )
    
    
    # the output table for the summary tab ------------------------------------
    
    output$mysummary <- DT::renderDataTable(
      my_summaries()[[1]],
      rownames = F,
      options = list(
        dom = "t",
        autoWidth = TRUE,
        columnDefs = list(list(
          width = '10%', targets = c(1, 1)
        ))
      )
    )
    
    # the output text for the summary tab -------------------------------------
    
    output$agg_info <- renderText(my_summaries()[[2]])
    
    # event reactive (reacts to the slider and to filters) for the OTU tab -------------------
    
    filt_OTU_table <- eventReactive(
      c(
        input$OTU,
        input$ssamplestable_rows_all,
        input$sample_agg_level,
        input$tax_agg_level
      ),
      {
        # transform in rel. ab.
        # provides a default when the input is not available
        # OTU_table_f <- NULL # data.frame(warning = "no table to display")
        # this works on the prevalence table generated in my summaries without
        # rarefaction
        if (!is.null(my_summaries()[[3]])) {
          OTU_table_f <-  my_summaries()[[3]] / rowSums(my_summaries()[[3]])
          top_taxa_to_keep <-
            ifelse(ncol(OTU_table_f) > input$OTU,
                   input$OTU,
                   ncol(OTU_table_f))
          top_taxa <-
            my_summaries()[[4]] %>%
            dplyr::slice(1:top_taxa_to_keep) %>%
            pull(label)
          OTU_table_f <-
            as.data.frame(OTU_table_f) %>%
            dplyr::select(all_of(top_taxa))
          other_column <-
            1 - rowSums(OTU_table_f)
          OTU_table_f <-
            OTU_table_f %>% rownames_to_column(var = "sample_name") %>%
            bind_cols(., "Other" = other_column) %>%
            mutate_if(is.numeric, round, digits = 4)
        }
        return(OTU_table_f)
      }
    )
    
    # diagnostic
    if (debug_mode) {
      output$sample_sums <- renderText(c(
        input$OTU,
        ncol(my_summaries()[[3]]) >
          input$OTU,
        rowSums(filt_OTU_table()[, 2:ncol(filt_OTU_table())])
      ))
    }
    
    if(debug_mode) {
      output$OTU_test <- DT::renderDataTable(FMBNexp()[[2]])
    }
    
    
    # the output table for the OTU tab ----------------------------------------
    
    output$OTU_table_agg <- DT::renderDataTable(
      filt_OTU_table(),
      rownames = F,
      options = list(dom = "lftip",
                     scrollX = T)
    )
    
    # the reactive for the save button in the OTU table tab -------------------
    
    observeEvent(input$save_OTU, {
      write_tsv(filt_OTU_table(),
                file.path(
                  "output",
                  "table_plots",
                  paste(
                    input$fn_prefix,
                    input$save_OTU,
                    "_filt_OTU_table.txt",
                    sep = ""
                  )
                ))
      showNotification(
        "filtered OTU table saved to /output/table_plots folder",
        type = "message",
        duration = 5,
        closeButton = TRUE
      )
    })
    
    # the output table for the prev & ab tab ----------------------------------
    
    output$prev_ab_df <- DT::renderDataTable(
      my_summaries()[[5]],
      rownames = F,
      options = list(dom = "lftip",
                     scrollX = T)
    )
    
    
    # the save button for the prev & ab table ---------------------------------
    
    observeEvent(input$save_prev_tab,
                 {
                   write_tsv(my_summaries()[[5]],
                             file.path(
                               "output",
                               "table_plots",
                               paste(
                                 input$fn_prefix,
                                 input$save_prev_tab,
                                 "_prev_ab_table.txt",
                                 sep = ""
                               )
                             ))
                   showNotification(
                     "filtered prevalence and abundance table saved to /output/table_plots folder",
                     type = "message",
                     duration = 5,
                     closeButton = TRUE
                   )
                 })
    
    
    # the prevalence and abundance plot ---------------------------------------
    
    # respond to both sliders and to changes in filter and aggregation
    tofilter <- eventReactive(
      c(
        input$ssamplestable_rows_all,
        input$sample_agg_level,
        input$tax_agg_level,
        input$prev,
        input$minab
      ),
      {
        # initialize with NULL
        # pass_prev_ab_filter <- NULL
        if (!is.null(my_summaries()[[5]])) {
          pass_prev_filter <- dplyr::filter(my_summaries()[[5]],
                                            Rel_Prev >= minprev) %>% pull(label)
          if (input$prev != minprev) {
            pass_prev_filter <- dplyr::filter(my_summaries()[[5]],
                                              Rel_Prev >= input$prev) %>% pull(label)
          }
          pass_prev_ab_filter <- dplyr::filter(my_summaries()[[5]],
                                               Max_Rel_Ab >= input$minab) %>% pull(label)
          # now the intersection
          pass_prev_ab_filter <-
            intersect(pass_prev_filter, pass_prev_ab_filter)
          
        }
        return(pass_prev_ab_filter)
      },
      ignoreNULL = F
    )
    
    
    # now the plot
    prev_ab_plot <- reactive({
      if (is.null(!is.null(my_summaries()[[5]]))) {
        NULL
      } else {
        prev_df <- my_summaries()[[5]]
        # first the default filter
        # prev_df <- prev_df %>%
        #   mutate(pass_prev_ab = ifelse((Rel_Prev>=input$prev && Max_Rel_Ab >= input$minab), "T", "F"))
        prev_df_f <- prev_df %>%
          mutate(pass_prev_ab = ifelse(label %in% tofilter(), "T", "F"))
        OTU0 <- base::ncol(my_summaries()[[3]])
        OTUmatrixf2 <- my_summaries()[[3]]
        OTUmatrixf2 <- OTUmatrixf2[, tofilter()]
        f_seq_ret <-
          round(sum(OTUmatrixf2) / sum(my_summaries()[[3]]), 4)
        title_text <- "Prevalence vs. abundance, by Phylum"
        taxa_counts <- prev_df_f %>%
          dplyr::filter(pass_prev_ab == "T") %>%
          summarise_at(.vars = c("phylum", "class", "order", "family", "genus"), n_distinct, na.rm=TRUE)
        subtitle_text <-
          str_wrap(paste(
            "Using the filters you retain ",
            length(tofilter()),
            " taxa (triangles) out of ",
            OTU0,
            " (",
            f_seq_ret * 100,
            "% of init. seqs.). ",
            "They belong to ",
            taxa_counts$phylum, " phyla, ",
            taxa_counts$class, " classes, ",
            taxa_counts$order, " orders, ",
            taxa_counts$family, " families, ",
            taxa_counts$genus, " genera.",
            sep = ""
          ), width = 80)
        # the plot
        prev_ab_plot <-
          ggplot(prev_df_f,
                 aes(
                   x = Tot_Ab,
                   y = Rel_Prev,
                   shape = as.factor(pass_prev_ab),
                   color = phylum
                 )) +
          geom_point(size = 2, alpha = 0.7) +
          facet_wrap(~ phylum) +
          geom_hline(
            yintercept = input$prev,
            alpha = 0.5,
            linetype = 2
          ) +
          labs(
            x = "total abundance",
            y = "Prevalence [Frac. Samples]",
            shape = 'pass ab. treshold',
            title = title_text,
            subtitle = subtitle_text
          ) +
          scale_x_log10() +
          scale_y_continuous(minor_breaks = seq(0, 1, 0.05)) +
          theme(
            legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90),
            strip.text = element_text(size = 6)
          )
        return(prev_ab_plot)
      }
    })
    
    
    output$prev_ab_plot <- renderPlot(if (is.null(my_summaries()[[5]])) {
      NULL
    } else {
      prev_ab_plot()
    })
    
    # diagnostic
    if (debug_mode) {
      output$slider_prev_value <- renderText(input$prev)
      output$slider_ab_value <- renderText(input$minab)
      output$pass_filter_text <- renderText(tofilter())
    }
    
    # the reactive for the save button in the prev. & ab. plot tab -------------------
    
    observeEvent(input$save_prev_graph, {
      ggsave(
        filename = file.path(
          "output",
          "table_plots",
          paste(
            input$fn_prefix,
            input$save_prev_graph,
            "_prev_ab_plot.",
            g_ext,
            sep = ""
          )
        ),
        
        plot = prev_ab_plot(),
        width = g_width,
        height = g_height,
        units = g_units,
        dpi = g_dpi
      )
      showNotification(
        "prev.&ab. plot saved to /output/table_plots folder",
        type = "message",
        duration = 5,
        closeButton = TRUE
      )
    })
    
    # the reactives for the bar plot ------------------------------------------
    
    # the dynamic menu for the sample aggregation level (bar and box plot tabs)
    
    sample_agg_menu_items <- reactive({
      # which categories in the sample table have more than one option?
      if (FMBNexp()[[11]] == "sample") {
        ncats <- FMBNexp()[[4]] %>%
          summarise_at(.vars =c("label", "foodId", "llabel", "L1", "L4", "L6"), n_distinct)
        c("none", "foodId", "exp. code", "L1", "L4", "L6")[which(ncats>1)]
      } else {
        ncats <- FMBNexp()[[4]] %>%
          summarise_at(.vars =c("label", "L1", "L4", "L6"), n_distinct)
        c("none", "L1", "L4", "L6")[which(ncats>1)]
      } 
    }
    ) # else is "exp. code"
    
    output$sample_agg_lvl <- renderUI(
      {
        selectInput(
          "sample_agg",
          label = "sample agg.",
          choices = sample_agg_menu_items(),
          multiple = FALSE
        )
      }
    )
    
    # diagnostic
    if (debug_mode) {
      output$sample_menu_items <- renderText(sample_agg_menu_items())
    }
    
    # the dynamic menu for the sample aggregation level (bar and box plot tabs)
    
    taxa_agg_menu_items <- reactive({
      if (input$tax_agg_level == "species") {
        c("none", colnames(FMBNexp()[[5]])[2:7])
      } else {
        c("none", colnames(FMBNexp()[[5]])[2:(which(colnames(FMBNexp()[[5]]) == FMBNexp()[[12]]) -
                                                1)])
      }
    } # else can be genus, family or class in the app)
    )
    
    output$taxa_agg_lvl <- renderUI({
      selectInput(
        "taxa_agg",
        label = "taxa agg.",
        choices = as.list(taxa_agg_menu_items()),
        selected = NULL,
        multiple = FALSE
      )
    })
    
    # diagnostic
    if (debug_mode) {
      output$taxa_menu_items <- renderText(taxa_agg_menu_items())
    }
    
    # the bar plot
    # the first reactive, generates the table used for the bar plot and box plot
    bb_plot_data <- eventReactive(
      c(
        input$ssamplestable_rows_all,
        input$sample_agg_level,
        input$sample_agg,
        input$prev,
        input$minab
      ),
      {
        # need to initialize with NULL?
        bb_plot_list_1 <- vector(mode = "list", length = 4)
        OTU_f_for_barplot <- my_summaries()[[3]]
        # calculate sums of sequences by sample for the "Other" column
        seq_sums_nofilter <- rowSums(my_summaries()[[3]])
        # apply the prev and ab filter
        OTU_f_for_barplot <- OTU_f_for_barplot[, tofilter()]
        n_removed <- ncol(my_summaries()[[3]])-ncol(OTU_f_for_barplot)
        # seq sums, filtered table
        seq_sums_filter <- rowSums(OTU_f_for_barplot)
        # difference
        Other_column <- seq_sums_nofilter-seq_sums_filter
        # add the column
        OTU_f_for_barplot <- cbind(OTU_f_for_barplot, Other = Other_column)
        # now melt and add the sample metadata for the first aggregation
        
        edge_table_filtered <- reshape2::melt(t(OTU_f_for_barplot))
        colnames(edge_table_filtered)[1:3] <- c("source", "target", "weight")
        
        if(length(sample_agg_menu_items()) == 4) {
          sel_s_var <- c("label", "L1", "L4", "L6")
        } else {
          sel_s_var <- c("label", "foodId", "llabel", "L1", "L4", "L6")
        }
        
        sample_metadata_j <- dplyr::select_at(FMBNexp()[[4]], .vars = sel_s_var) %>%
          dplyr::rename(target = label)
        
        edge_table_filtered_j <- left_join(edge_table_filtered, sample_metadata_j)
        
        bxp_edge_table <- edge_table_filtered_j
        
        sum_abs <- bxp_edge_table %>%
          group_by(target) %>%
          dplyr::summarise(sumab = sum(weight, na.rm = T))
        
        bxp_edge_table_r <- left_join(bxp_edge_table, sum_abs) %>%
          mutate(relab = weight/sumab) %>%
          mutate(relab_c = relab+0.00001)
        
        # translating the select input; the name is the menu choice, the value is the actual variable
        sample_menu_choice <- input$sample_agg
        
        if(sample_menu_choice  == "none") {
          if(FMBNexp()[[11]] == "sample"){
            sample_agg_choice <- "label"
            names(sample_agg_choice) <- FMBNexp()[[11]]
          } else {
            sample_agg_choice <- "llabel"
            names(sample_agg_choice) <- FMBNexp()[[11]]  
          }
        } else {
          if(sample_menu_choice  == "exp. code") {
            sample_agg_choice <- "llabel"
            names(sample_agg_choice) <- "exp. code"
          } else {        
            sample_agg_choice <- sample_menu_choice
            names(sample_agg_choice) <- sample_menu_choice}
        }
        # that is: if "none" the variable used for joining is "label" but the name is either "sample" or "exp. code"
        # now need to perform a group_by_at
        # determine which of the sample group choices result in a number of items <= max_samples (see options)
        sgroup_items <- sort(sapply(FMBNexp()[[4]][, sel_s_var], n_distinct, na.rm = T))
        sgroup_items_pass <- names(sgroup_items)[which(sgroup_items <= max_samples)]
        
        sample_grouping <- ifelse(
          sample_agg_choice %in% sgroup_items_pass,
          sample_agg_choice, 
          sgroup_items_pass[length(sgroup_items_pass)]
        )
        # returns "label" when aggr is exp. code and number of llabels <25
        
        sgrouping_message_flag <- F # initialize the grouping message flag
        
        # do not aggregate if FMBNexp()[[13]] == T
        if(names(sample_agg_choice) == FMBNexp()[[11]] | FMBNexp()[[13]]) {
          if ("label" %in% sgroup_items_pass) {
            edge_table_filtered_j_s  <- edge_table_filtered_j[, 1:3] %>%
              dplyr::rename(sweight = weight)
          } else {
            edge_table_filtered_j_s  <- edge_table_filtered_j %>%
              group_by_at(.vars = c("source", sample_grouping)) %>%
              summarise(sweight = sum(weight, na.rm = T))
            sgrouping_message_flag <- T
            sgrouping_message <-
              paste("WARNING: too many categories, grouping by ",
                    sample_grouping,
                    " instead!",
                    sep = "")
            
            showNotification(
              sgrouping_message,
              type = "message",
              duration = 10,
              closeButton = TRUE
            )
          }
        } else {
          if (sample_agg_choice %in% sgroup_items_pass) {
            edge_table_filtered_j_s  <- edge_table_filtered_j %>%
              group_by_at(.vars = c("source", sample_agg_choice)) %>%
              summarise(sweight = sum(weight, na.rm = T))
          } else {
            edge_table_filtered_j_s  <- edge_table_filtered_j %>%
              group_by_at(.vars = c("source", sample_grouping)) %>%
              summarise(sweight = sum(weight, na.rm = T))
            sgrouping_message_flag <- T
            sgrouping_message <-
              paste("WARNING: too many categories, grouping by ",
                    sample_grouping,
                    " instead!",
                    sep = "")
            showNotification(
              sgrouping_message,
              type = "message",
              duration = 10,
              closeButton = TRUE
            )
          }
        }
        
        colnames(edge_table_filtered_j_s)[2] <- "target"
        taxa_sel <- c("label", as.character(taxa_agg_menu_items()[-1]))
        taxa_j <- dplyr::select_at(FMBNexp()[[5]], .vars = taxa_sel) %>%
          dplyr::rename(source = label)
        # join taxa info
        edge_table_filtered_j_s_t <- left_join(edge_table_filtered_j_s, taxa_j) 
        
        bxp_edge_table_r_j <- left_join(bxp_edge_table_r, taxa_j)
        
        # for diagnostic purposes the first element can be shown as a table
        # change the first element of the list:
        # alternatives are 
        # my_summaries()[[3]] the OTU table from my summaries
        # FMBNexp()[[4]] the sample metadata table
        # OTU_f_for_barplot the filtered OTU table (with the Other column)
        # edge_table_filtered the melted table 
        # sample_metadata_j the sample metadata
        # edge_table_filterd_j the edge table with metadata
        # edge_table_filtered_j_s the edge table with metadata after aggregation by sample_agg
        # FMBNexp()[[5]] the taxa metadata after aggregation (including NA columns)
        # taxa_j the taxa metadata, with label renamed and just the colums in the tax agg menu
        # edge_table_filtered_j_s_t edge table + taxa info (used by the box plot)
        # bxp_edge_table_r the object for the box plot
        bb_plot_list_1[[1]] <- bxp_edge_table_r_j
        bb_plot_list_1[[2]] <- edge_table_filtered_j_s_t
        bb_plot_list_1[[3]] <- bxp_edge_table_r_j
        bb_plot_list_1[[4]] <- sample_grouping
        bb_plot_list_1[[5]] <- taxa_j
        return(bb_plot_list_1)
      },
      ignoreNULL = F
    )
    
    # the reactive for the bar plot, responds to the taxa aggregation menu but also to
    # changes in other menus
    bb_plot_data_2 <- eventReactive(
      c(
        input$ssamplestable_rows_all,
        input$tax_agg_level,
        input$sample_agg_level,
        input$sample_agg,
        input$taxa_agg,
        input$prev,
        input$minab
      ),{
        
        bb_plot_list_2 <- vector(mode = "list", length = 3)
        input_egde_table <- bb_plot_data()[[2]]
        # determine which of the taxa group choices result in a number of items <= max_taxa (see options)
        sel_t_vars <- c("source", as.character(taxa_agg_menu_items())[-1])
        tgroup_items <- sapply(input_egde_table[, sel_t_vars], n_distinct, na.rm = T)
        tgroup_items_pass <- names(tgroup_items)[which(tgroup_items<=max_taxa)]
        taxa_grouping <- ifelse(
          input$taxa_agg %in% tgroup_items_pass,
          input$taxa_agg, 
          tgroup_items_pass[length(tgroup_items_pass)]
        )
        
        tgrouping_message_flag <- F # initialize the grouping message flag
        
        if(input$taxa_agg == "none") {
          if ("source" %in% tgroup_items_pass) {
            edge_table  <- input_egde_table[, 1:3] %>%
              dplyr::rename(weight = sweight)
          } else {
            edge_table  <- input_egde_table %>%
              group_by_at(.vars = c(taxa_grouping, "target")) %>%
              summarise(weight = sum(sweight, na.rm = T)) %>%
              ungroup()
            tgrouping_message_flag <- T
            tgrouping_message <-
              paste("WARNING: too many taxa, grouping by ",
                    taxa_grouping,
                    " instead!",
                    sep = "")
            showNotification(
              tgrouping_message,
              type = "message",
              duration = 10,
              closeButton = TRUE
            )
          }
        } else {
          if (input$taxa_agg %in% tgroup_items_pass) {
            edge_table  <- input_egde_table %>%
              group_by_at(.vars = c(input$taxa_agg, "target")) %>%
              summarise(weight = sum(sweight, na.rm = T)) %>%
              ungroup()
          } else{
            edge_table  <- input_egde_table %>%
              group_by_at(.vars = c(taxa_grouping, "target")) %>%
              summarise(weight = sum(sweight, na.rm = T)) %>%
              ungroup()
            tgrouping_message_flag <- T
            tgrouping_message <-
              paste("WARNING: too many taxa, grouping by ",
                    taxa_grouping,
                    " instead!",
                    sep = "")
            showNotification(
              tgrouping_message,
              type = "message",
              duration = 10,
              closeButton = TRUE
            )
          }
        }
        
        # number of distinct taxa
        n_taxa <- n_distinct(edge_table[, 1], na.rm = F)
        names(edge_table) <- c("source", "target", "abs_ab")
        # transform NA in Other and sum if necessary, recode as a factor and move "Other"
        # at the end if necessary, to provide a nicer graph
        if(anyNA(edge_table$source)){
          edge_table <- edge_table %>%
            mutate(source = if_else(is.na(source),"Other",source)) %>%
            group_by(source, target) %>%
            summarise(abs_ab = sum(abs_ab, na.rm = T)) %>%
            ungroup() %>%
            mutate(source = factor(source, ordered = T)) 
          if("Other" %in% levels(edge_table$source)){
            edge_table <- edge_table %>% 
              mutate(source = forcats::fct_relevel(source, "Other", after = 0))
          }
        }
        
        # the bar plot
        bar_plot_title <- paste("Relative abundance,", 
                                taxa_grouping, 
                                sep = " ")
        bar_plot_subtitle <- str_wrap(
          paste("Only taxa with >",
                input$prev,
                "prevalence and >",
                input$minab,
                "abundance are shown. The others are aggregated as \"Other\".",
                sep = " "),
          width = 80
        )
        if (n_taxa>12) {
          gpalette <- rpalette[1:n_taxa] 
          bar_plot <- ggplot(edge_table, mapping = aes(x = target, y = abs_ab, fill = source)) + 
            geom_col(position = "fill") + 
            scale_fill_manual(values = gpalette) +
            labs(title = bar_plot_title,
                 subtitle = bar_plot_subtitle,
                 x= "Samples/Sample groups",
                 y= "Relative abundance",
                 fill = "Taxa") +
            scale_x_discrete(labels = function(x) str_wrap(x, str_wrap_length)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                  panel.background= element_rect(fill = "white"),
                  plot.title = element_text(hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5, size = 9))
        } else {
          bar_plot <- ggplot(edge_table, mapping = aes(x = target, y = abs_ab, fill = source)) + 
            geom_col(position = "fill") + 
            labs(title = bar_plot_title,
                 subtitle = bar_plot_subtitle,
                 x= "Samples/Sample groups",
                 y= "Relative abundance",
                 fill = "Taxa") +
            scale_fill_brewer(type = "qual", palette = "Paired") +
            scale_x_discrete(labels = function(x) str_wrap(x, str_wrap_length)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                  panel.background= element_rect(fill = "white"),
                  plot.title = element_text(hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5, size = 9))
        }
        
        bb_plot_list_2[[1]] <- edge_table # also used as diagnostic
        bb_plot_list_2[[2]] <- n_taxa
        bb_plot_list_2[[3]] <- bar_plot
        bb_plot_list_2[[4]] <- tgroup_items_pass
        return(bb_plot_list_2)  
      },
      ignoreNULL = F
      
    )
    
    # for diagnostic use
    output$barplot_table <- if(debug_mode) {DT::renderDataTable(bb_plot_data()[[1]])} # or _1
    
    output$barplot <- renderPlot(bb_plot_data_2()[[3]])
    
    
    # the reactive for saving the bar plot
    observeEvent(input$save_bar_plot, {
      ggsave(
        filename = file.path(
          "output",
          "table_plots",
          paste(
            input$fn_prefix,
            as.integer(input$save_bar_plot),
            "_bar_plot.",
            g_ext,
            sep = ""
          )
        ),
        
        plot = bb_plot_data_2()[[3]],
        width = g_width,
        height = g_height,
        units = g_units,
        dpi = g_dpi
      )
      showNotification(
        "bar plot saved to /output/table_plots folder",
        type = "message",
        duration = 5,
        closeButton = TRUE
      )
    })
    
    
    # the reactives for the box plot ------------------------------------------
    
    # the dynamic menu for the taxa categories
    output$taxa_agg_bxp <- renderUI({
      selectInput(
        "taxa_agg_bx",
        label = "taxa agg.",
        choices = as.list(bb_plot_data_2()[[4]][-1]),
        selected = "2",
        multiple = FALSE
      )
    })
    
    # the reactive for the taxa item menu
    taxa_items <- eventReactive(
      c(
        input$ssamplestable_rows_all,
        input$tax_agg_level,
        input$sample_agg_level,
        input$sample_agg,
        input$taxa_agg,
        input$prev,
        input$minab,
        input$taxa_agg_bx,
        input$g_type,
        input$taxa_agg_bx
      ),{
        my_vars <- c("source", input$taxa_agg_bx)
        my_items <- select_at(bb_plot_data()[[2]], .vars = my_vars)
        my_item_list <- dplyr::distinct_at(my_items, .vars = input$taxa_agg_bx) %>%
          select_at(.vars = input$taxa_agg_bx) %>%
          pull()
        return(my_item_list[complete.cases(my_item_list)])
      }
    )
    
    # the contextual menu with the taxa items
    output$sel_taxa <- renderUI({
      selectInput(
        "sel_taxa_items",
        label = "select taxa",
        choices = c("all",taxa_items()),
        selected = "all",
        multiple = TRUE
      )
    })
    
    # the list containing data and boxplot
    # NOTE it might make sense to add an option for selecting an arithmetic scale
    bb_plot_data_3 <- eventReactive(
      c(
        input$ssamplestable_rows_all,
        input$tax_agg_level,
        input$sample_agg_level,
        input$sample_agg,
        input$taxa_agg,
        input$prev,
        input$minab,
        input$taxa_agg_bx,
        input$sel_taxa_items,
        input$g_type
      ),
      {
        bb_plot_list_3 <- vector(mode = "list", length = 3)
        
        bxp_starting_et <- bb_plot_data()[[3]]
        
        # the sample aggregation level (max 25 cats. and can't be sample), this handles the case
        # when no aggregation is possible. Depends on the menu in bar plot table
        sample_agg_bp <- ifelse(bb_plot_data()[[4]] == "label", "target", bb_plot_data()[[4]])
        
        
        # the taxa aggregation options (domain excluded)
        taxa_item_names <- bb_plot_data_2()[[4]][-1] 
        
        # selected taxa category
        sel_taxa_cat <- ifelse(is.null(input$taxa_agg_bx),"class",input$taxa_agg_bx)
        
        # select variables and remove NAs
        bxp_starting_et_2 <- bxp_starting_et %>%
          select_at(.vars = c(sel_taxa_cat, sample_agg_bp, "relab", "relab_c")) 
        colnames(bxp_starting_et_2) <- c("taxa", "samples", "relab", "relab_c")
        bxp_starting_et_2 <- dplyr::filter(bxp_starting_et_2, !is.na(taxa))
        
        # now filter using the input$sel_taxa_items
        if("all" %in% input$sel_taxa_items){
          bxp_starting_et_3 <- bxp_starting_et_2 
        } else {
          t_selection <- input$sel_taxa_items[which(input$sel_taxa_items != "all")]
          if(length(t_selection)<2){
            bxp_starting_et_3 <- bxp_starting_et_2
            showNotification(
              "you must select at least two taxa; all will be selected, try again",
              type = "message",
              duration = 5,
              closeButton = TRUE
            )
          } else {
            bxp_starting_et_3 <- bxp_starting_et_2 %>%
              dplyr::filter(taxa %in% t_selection)
          }
        }
        
        
        # the actual box plot, first by selecting everything given the category menu
        box_plot_title <- paste("Abundance distribution,", 
                                sel_taxa_cat, 
                                sep = " ")
        
        
        # handle rare case when sample_agg_bp == "label"
        if(sample_agg_bp == "target"){
          top_x_boxplot <- ggplot(bxp_starting_et_3, mapping = aes(x = taxa, y = relab_c))
          if(input$g_type == "violin"){
            top_x_boxplot <- top_x_boxplot +
              geom_violin() +
              geom_jitter(alpha = 0.2, width = 0.2) +
              labs(x= "Taxa", 
                   y= "rel.ab.") +
              scale_x_discrete(labels = function(x) str_wrap(x, 15)) +
              scale_y_log10(breaks = c(1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1)) +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                    strip.text = element_text(size = 6)) 
          } else {
            top_x_boxplot <- top_x_boxplot +
              geom_boxplot() +
              geom_jitter(alpha = 0.2, width = 0.2) +
              labs(x= "Taxa", 
                   y= "rel.ab.") +
              scale_x_discrete(labels = function(x) str_wrap(x, 15)) +
              scale_y_log10(breaks = c(1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1)) +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                    strip.text = element_text(size = 6)) 
          }
        } else {
          top_x_boxplot <- ggplot(bxp_starting_et_3, mapping = aes(x = samples, y = relab_c))
          if(input$g_type == "violin"){
            top_x_boxplot <- top_x_boxplot +
              geom_violin() +
              geom_jitter(alpha = 0.2, width = 0.2) +
              facet_wrap(~taxa) +
              labs(x= "Samples/Sample groups", 
                   y= "rel.ab.") +
              scale_x_discrete(labels = function(x) str_wrap(x, 15)) +
              scale_y_log10(breaks = c(1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1)) +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                    strip.text = element_text(size = 6)) 
          } else {
            top_x_boxplot <- top_x_boxplot +
              geom_boxplot() +
              geom_jitter(alpha = 0.2, width = 0.2) +
              facet_wrap(~taxa) +
              labs(x= "Samples/Sample groups", 
                   y= "rel.ab.") +
              scale_x_discrete(labels = function(x) str_wrap(x, 15)) +
              scale_y_log10(breaks = c(1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1)) +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                    strip.text = element_text(size = 6)) 
          }   
        }
        
        
        # make the list
        bb_plot_list_3[[1]] <- bxp_starting_et_3 # for diagnostic purposes
        bb_plot_list_3[[2]] <- top_x_boxplot
        
        return(bb_plot_list_3)  
        
      }
    )
    
    # diagnostic for box plot
    output$boxplot_table <- if(debug_mode) {DT::renderDataTable(bb_plot_data_3()[[1]])} 
    # the boxplot
    output$topxboxplot <- renderPlot(bb_plot_data_3()[[2]])
    
    # the reactive for saving the box plot
    observeEvent(input$save_box_plot, {
      ggsave(
        filename = file.path(
          "output",
          "table_plots",
          paste(
            input$fn_prefix,
            as.integer(input$save_box_plot),
            "_box_plot.",
            g_ext,
            sep = ""
          )
        ),
        
        plot = bb_plot_data_3()[[2]],
        width = g_width,
        height = g_height,
        units = g_units,
        dpi = g_dpi
      )
      showNotification(
        "box plot saved to /output/table_plots folder",
        type = "message",
        duration = 5,
        closeButton = TRUE
      )
    })
      
  #   closing curly bracket for server    
    
  } 
  
  
  
  ## this is for loading helpers, if any. They must be in the same folder of the app
  # source("helpers.R")
  
  # Create Shiny app ----
  shinyApp(ui = ui, server = server)
