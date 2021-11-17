# ShinyFMBN v1.2-----------------------------------------------------

# a Shiny app to explore, filter and extract data from FoodMicrobionet

# preparatory steps -------------------------------------------------------

  # install/load packages -------------------------------------------------

  .cran_packages <- c("shiny", "DT", "tidyverse", "reshape2", "stringr", 
                      "tibble", "lubridate" ,"igraph")
  .bioc_packages <- c("BiocManager", "phyloseq")
  .inst <- .cran_packages %in% installed.packages()
  if(any(!.inst)) {
    install.packages(.cran_packages[!.inst])
  }
  .inst <- .bioc_packages %in% installed.packages()
  if(any(!.inst)) {
    if(!.inst[1]) install.packages("BiocManager")
    if(any(!.inst[2:length(.inst)])) {
      BiocManager::install(!.inst[2:length(.inst)], ask = F)
    }
  }
  sapply(c(.cran_packages, .bioc_packages), require, 
         character.only = TRUE) 
  
  # load files (from the data folder) -------------------------------------
  # using a platform independent path
  studies <- readRDS(file.path("data","studies.rds"))
  samples <- readRDS(file.path("data","samples.rds"))
  edges <- readRDS(file.path("data", "edges.rds"))
  taxa <- readRDS(file.path("data", "taxa.rds"))
  references <- readRDS(file.path("data", "references.rds"))
  version <- readRDS(file.path("data", "version.rds"))
  
  # generate minitables for display
  studies_disp <- studies %>% dplyr::select(studyId, FMBN_version, target, region, 
                                     platform, Seq_accn, Seq_accn_link, 
                                     short_descr, samples, ref_short, DOI, DOI_link)
  studies_disp <- studies_disp %>% 
    mutate(
      target = paste(target, ", ", region, sep =""),
      SRA_ENA = ifelse(is.na(Seq_accn),
                       "not available", 
                       paste0("<a href='",
                              Seq_accn_link,
                              "' target='_blank'>",
                              Seq_accn,"</a>")),
      DOIlink = ifelse(DOI == "unpublished data",
                       "unpublished data",
                       paste0("<a href='",
                              DOI_link,
                              "' target='_blank'>",
                              DOI,"</a>"))
    ) %>%
    dplyr::select(studyId, FMBN_version, target, platform, samples, SRA_ENA,
           ref_short, DOIlink, short_descr)
  
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
  
  FMBN_summary_text <- paste("There are", n_studies, "studies and", n_samples, 
                             "samples in this version of FoodMicrobionet.", 
                             "The samples belong to", n_food_groups,
                             "major food groups and", n_foodId, "different foods.",
                             "There are", n_llabel, 
                             "different combinations of food, nature, process, fermentation/spoilage.",
                             n_taxa, "taxa have been identified at different taxonomic levels.")
  
# user interface ----------------------------------------------------------

ui <- navbarPage("Shiny FMBN",
                 # View tab --------------------------------------------------------------
                 tabPanel("Explore",
                          fluidPage(
                            # row for studies table
                            fluidRow(
                              column(3, 
                                     h4("FMBN studies"),
                                     p("Explore studies in FMBN by navigating the table on the right."),
                                     tags$div(
                                       tags$ul(
                                         tags$li(tags$span("Use the search box to find studies by keywords in any field")),
                                         tags$li(tags$span(
                                           "Use the links to access studies on SRA/ENA or articles via DOI."
                                         )),
                                         tags$li(tags$span("Hover with the mouse over a description cell to get a full description."))
                                         )
                                       )),
                              column(9, div(DT::dataTableOutput("studies_table"), style = "font-size:80%"))
                            ),
                            # row for summary stats on FMBN
                            fluidRow(
                              column(3, 
                                     hr(),
                                     h4("Summary statistics on FMBN")),
                              column(9, 
                                     hr(),
                                     p(paste(version, FMBN_summary_text)))
                            ),
                            # row for selecting and viewing a study
                            fluidRow(
                              column(3, 
                                     hr(),
                                     h4("Samples"),
                                     uiOutput('choose_study')
                               ),
                                     
                              column(9, 
                                     hr(),
                                     div(DT::dataTableOutput("vsamples_table"), style = "font-size:80%")
                                     )
                              
                            )
                          )),
                 # Filter tab --------------------------------------------------
                 tabPanel("Filter",
                          fluidPage(
                  # Filter instructions  and widgets ----------------------------          
                            fluidRow(
                              column(2,
                                     selectInput("select_study", label = "Select a subset of studies", 
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
                                     p("Use dropdown menus to select by study, food group, food code (FoodEx2). 
                                       Selected samples will appear in the table, which can be used  
                                       as a guidance for refining your selection.",
                                       "Dropdown filters are applied in the order they appear in this page. 
                                       Use the <apply> checkbox to apply/remove a filter from a drop down menu. ",
                                       "When you are done, you can use the filters at the bottom of the table
                                       (sample type, nature, process, target) to further refine your search.", 
                                       style = "font-size:80%")
                                      ),
                  # Filtered table ----------------------------------------------
                              column(10, 
                                     div(DT::dataTableOutput("ssamplestable"), style = "font-size:80%")
                                     )
                            )

                          )),
                 # Aggregate tab -----------------------------------------------------------
                 tabPanel("Aggregate",
                          fluidPage(
                            fluidRow(
                              column(2,
                                     selectInput("sample_agg_level", label = "Select the sample aggregation level", 
                                                 choices = c("sample", "exp. code"), 
                                                 selected = "sample",
                                                 multiple = F
                                                 ),
                                     selectInput("tax_agg_level", label = "Select the taxonomic aggregation level", 
                                                 choices = c("species", "genus", "family", "class"), 
                                                 selected = "species",
                                                 multiple = F
                                                 ),
                                     selectInput("show_what", label = "Show", 
                                                 choices = c("original", "aggregated"), 
                                                 selected = "original",
                                                 multiple = F
                                     ),
                                     hr(),
                                     strong("Instructions"),
                                     p("Once you have selected a subset of samples in the Filter tab, ",
                                       "tables for studies, samples, taxa (in descending order of ",
                                       "relative abundance), and references will appear on the left.", style = "font-size:80%"),
                                     p("If you are happy with your selection, use the tools in this tab ",
                                       "to (optionally) aggregate samples and taxa and to show original or ",
                                       "aggregated data and then move to the Export tab to save your files. ",
                                       "Otherwise go back to the Filter tab and start over.", style = "font-size:80%"),
                                     strong("Please note:"),
                                     p("When combining studies using different platforms/pipelines, you should choose",
                                       "genus or above as taxonomic aggregation level.", style = "font-size:80%"),
                                     p("When you choose the expanded food code as a sample aggregation level, ",
                                       "consider the risks related to aggregating samples for which different",
                                       "wet and dry lab procedure and bionformatics pipelines were used. ",
                                       "In addition at the latter aggregation level a heuristic is used to ",
                                       "recalculate the number of sequences for each group.", style = "font-size:80%"),
                                     p("Refreshing the display for large tables may take some time.", style = "font-size:80%")
                                     ),
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
                              )
                            )
                          ),
                 # Export tab --------------------------------------------------------------
                 tabPanel("Export",
                          fluidPage(
                            fluidRow(
                              column(3,
                                     textInput("fn_prefix", "Choose a filename prefix", value = "myFMBN"),
                                     hr(),
                                     h5("Export your files."),
                                     actionButton("miniFMBN","urFMBN"),
                                     actionButton("aggFMBN","aggFMBN"),
                                     actionButton("phyloseq","phyloseq"),
                                     actionButton("conet","conet"),
                                     actionButton("gml","gml"),
                                     hr(),
                                     strong("Instructions"),
                                     p("Summary data on your selection, a scrollable OTU abundance table ",
                                       " from your selection and a status message will appear in the right pane.", 
                                       style = "font-size:80%"),
                                     p("Use the text input to provide a prefix which will be included ",
                                       "in all your filenames.", style = "font-size:80%"),
                                     p("Use the action buttons to export files extracted from FoodMicrobionet in the output folder:",
                                       style = "font-size:80%"),
                                     tags$div(
                                       tags$ul(
                                         tags$li(tags$span(
                                           "urFMBN: a .rds file with filtered FMBN tables, no aggregation (subfolder minifmbn)"
                                         )),
                                         tags$li(tags$span(
                                           "aggFMBN: a .rds file with FMBN tables after aggregation, for further graphical and statistical analysis (subfolder aggdata)"
                                         )),
                                         tags$li(tags$span("phyloseq: a phyloseq class object, with OTU, sample and taxa tables, after aggregation, for use with Shiny-Phyloseq (subfolder phyloseq)")),
                                         tags$li(tags$span("conet: OTU and sample tables, after agregation, in a format suitable for use with the CoNet app (subfolder conet)")),
                                         tags$li(tags$span(
                                           "gml: the bipartite sample and taxa network in .gml format, can be imported in Cytoscape and Gephi (subfolder gml)"
                                         ))
                                       ),
                                       style = "font-size:80%"),
                                     p("A notification will briefly appear at the bottom right every time you export a file.",
                                       style = "font-size:80%"),
                                     p("Refreshing the display for large tables may take some time.", style = "font-size:80%")
                                     ),
                              column(9,
                                     h4("Summary info on your selection"),
                                     div(DT::dataTableOutput("mysummary"), style = "font-size:80%"),
                                     hr(),
                                     textOutput("agg_info"),
                                     hr(),
                                     h4("Your OTU table."),
                                     div(DT::dataTableOutput("OTU_table_agg"), style = "font-size:80%")
                                     )
                              
                            )
                            )
                          ),
                 # About tab ---------------------------------------------------------------
                 tabPanel("About",
                          fluidPage(
                            titlePanel("About"),
                            sidebarLayout(
                              sidebarPanel(
                                h4("About FoodMicrobionet"),
                                p("FoodMicrobionet is a database of studies on the
                                  composition of the bacterial microbiota of food obtained by
                                  amplicon targeted (16S RNA gene or 16S RNA) high throughput
                                  sequencing."
                                ),
                                p("Data are organized in tables:"),
                                tags$div(
                                  tags$ul(
                                    tags$li(tags$span("Studies: metadata on studies")),
                                    tags$li(tags$span(
                                      "Samples: metadata on samples for each study"
                                    )),
                                    tags$li(tags$span("Edges: taxa abundance for each sample")),
                                    tags$li(tags$span("Taxa: metadata on taxa"))
                                    )
                                  ),
                                p("To learn more on FoodMicrobionet and its development visit ",
                                  span(a("our homepage. ",
                                         href = "http://www.foodmicrobionet.org",
                                         target = "_blank")
                                       ),
                                  "To cite FoodMicrobionet and see related articles ",
                                  span(a("go here.",
                                         href = "https://doi.org/10.1016/j.ijfoodmicro.2015.12.001",
                                         target = "_blank")
                                       )
                                  ),
                                h4("About this app"),
                                p("ShinyFMBN1_2 is designed to provide an interface to:"),
                                tags$div(tags$ul(
                                  tags$li(tags$span("explore studies and samples;")),
                                  tags$li(tags$span("extract/filter groups of samples;")),
                                  tags$li(
                                    tags$span("(optionally) perform aggregation of samples or taxa;")
                                  ),
                                  tags$li(tags$span("(locally) save data for further analysis"))
                                  )
                                  ),
                                h4("Credits and copyright"),
                                p("The function for installing/loading packages is taken verbatim from ",
                                  span(
                                    a("https://f1000research.com/articles/5-1492/v2",
                                      href = "https://f1000research.com/articles/5-1492/v2",
                                      target = "_blank")
                                    )
                                  ),
                                h5("Copyright 2018 Eugenio Parente, Università della Basilicata"),
                                p("Permission is hereby granted, free of charge, to any person obtaining a copy ",
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
                              mainPanel(
                                p("A legacy image from FoodMicrobionet 2.0, generated with ",
                                  span(a("Gephi.", href = "https://gephi.org",
                                         target = "_blank")
                                       )
                                  ),
                                img(
                                  src = "FMBN2_0_small.jpg",
                                  height = 600,
                                  width = 600
                                  )
                                )
                              ) # closing bracket sidebar layout
                            ) # closing bracket fluidPage UI function
                          )
) # closing bracket UI function


# define server -----------------------------------------------------------

server <- function(input, output, session) {
 
  # non reactive outputs -----------------------------------------------------  
    # View tab, studies -------------------------------------------------------
    output$studies_table <- DT::renderDataTable(studies_disp, 
                                             rownames = F,
                                             escape = c(7,9),
                                             colnames = c("study" = "studyId", 
                                                          "version" = "FMBN_version", 
                                                          "SRA/ENA" = "SRA_ENA",
                                                          "reference" = "ref_short", 
                                                          "DOI link" = "DOIlink",
                                                          "description" = "short_descr"),
                                             options = list(pageLength = 3,
                                                            lengthMenu = c(3, 5, 10, 15),
                                                            columnDefs = list(
                                                              list(
                                                                targets = 8,
                                                                render = JS(
                                                                  "function(data, type, row, meta) {",
                                                                  "return type === 'display' && data.length > 100 ?",
                                                                  "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;","}"
                                                                  )
                                                                )
                                                              )
                                                            )
                                             )
  
  # samples table, view tab -------------------------------------------------
  # Render selectInput for study, View tab 
    
    output$choose_study <- renderUI({
    study_names <- studies$studyId
    selectInput("view_study",
                label = h5("Select one or more studies to explore samples"), 
                choices = study_names, 
                selected = "ST1",
                multiple = TRUE)    
  })
  
  # Reactively subset so that only the selected study is viewed

  f_vsamples <- reactive({
    subset(samples, studyId %in% input$view_study, 
           select = c(label_2, s_type, n_reads2, foodId, llabel, L1, description))
  })

  output$vsamples_table <- DT::renderDataTable(
    f_vsamples(),
    rownames = F,
    escape = T,
    colnames = c("label" = "label_2", 
                 "sample type" = "s_type", 
                 "reads" = "n_reads2",
                 "code" = "foodId", 
                 "ext. code" = "llabel",
                 "food group" = "L1",
                 "description" = "description"),
    options = list(pageLength = 5,
                   lengthMenu = c(5, 10, 25, 20),
                   columnDefs = list(
                     # limits display length of column 5 = description
                     list(
                       targets = 5,
                       render = JS(
                         "function(data, type, row, meta) {",
                         "return type === 'display' && data.length > 30 ?",
                         "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;","}"
                       )
                     )
                   )
                   
    )
  )
  
  # functions for the filter tab --------------------------------------------

  # gets the list of studies to be selected
  study_sel <- reactive({ 
                 if(is.null(input$select_study) | input$study_ch_box == F){
                   studies$studyId
                   } else {
                     input$select_study
                     }
                },
                 label = "study_sel_reactive")
  
  # reactive for the dynamic food group list
  food_group_names <- reactive(
    {
      nfstudy <- subset(samples, studyId %in% study_sel())
      unique(nfstudy$L1)
    }
  )
  
  # renders the food group drop down 
  output$choose_food_group <- renderUI({
    selectInput("select_L1",
                label = "Select one or more food groups", 
                choices = food_group_names(), 
                selected = NULL,
                multiple = TRUE)    
  })

  # creates the reactive which responds to select_L1
  food_group_sel <- reactive({ 
    if(is.null(input$select_L1) | input$food_group_ch_box == F){
      # using dplyr::pull to pull a vector
      food_group_names()
    } else {
      input$select_L1
    }
  },
  label = "food_group_sel_reactive")
  
  # reactive, for the dynamic foodId list
  
  foodId_names <- reactive(
    {
      nfstudy <- subset(samples, studyId %in% study_sel())
      nffgroup <- subset(nfstudy, L1 %in% food_group_sel())
      unique(nffgroup$foodId)
    }
  )
 
  # renders the foodId drop down
  output$choose_foodId <- renderUI({
    selectInput("select_foodId",
                label = "Select one or more food codes", 
                choices = foodId_names(), 
                selected = NULL,
                multiple = TRUE)    
  }) 

  # creates the reactive which responds to select_foodId
  foodId_sel <- reactive({ 
    if(is.null(input$select_foodId) | input$foodId_ch_box == F){
      foodId_names()
    } else {
      input$select_foodId
    }
  },
  label = "foodId_sel_reactive")

  # creates the reactive to display the table
  f_fsamples <- reactive(
    {
      fstudy <- subset(ssamples, studyId %in% study_sel())
      ffgroup <- subset(fstudy, L1 %in% food_group_sel())
      ffoodId <- subset(ffgroup, foodId %in% foodId_sel())
      ffinal <- ffoodId
      ffinal <- mutate(ffinal, SRA_run = ifelse(is.na(SRA_run),
                       "not available", 
                       paste0("<a href='",
                              "https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=",
                              SRA_run,
                              "' target='_blank'>",
                              SRA_run,"</a>")))
      return(ffinal)
    }
  )
  
  # renders the table for the filter tab 
  output$ssamplestable <- DT::renderDataTable(
    f_fsamples(),
    filter = "bottom",
    rownames = F,
    escape = 21,
    colnames = c("study" = "studyId",
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
                 "run" = "SRA_run"),
    options = list(pageLength = 10,
                   lengthMenu = c(10, 20, 30, 40),
                   scrollX = T,
                   columnDefs = list(
                      # limits the variables which are visible
                     list(visible=FALSE, targets=c(1,2,4,5,7,13,20,21)),
                     # defines which columns are searchable
                     list(targets = c(0,3,10,11,12,14,22), searchable = FALSE),
                      # limits number of characters shown, more on hover
                      list(targets = c(12,13,15),
                      render = JS(
                        "function(data, type, row, meta) {",
                        "return type === 'display' && data.length > 20 ?",
                        "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;","}"
                      ))
                   )
                   
    )
  )
  
  # functions for the aggregate tab -----------------------------------------
  
  FMBNfilt <- eventReactive(
    input$ssamplestable_rows_all, {
      # extract samples based on row numbers of the filtered table
      samples_filt <- samples %>% 
        dplyr::slice(
          pull(f_fsamples()[input$ssamplestable_rows_all,2])
        )
      # get studies
      studies_filt <- studies %>% 
        dplyr::filter(studyId %in% pull(distinct(dplyr::select(samples_filt, studyId), studyId)))
      # get references
      reference_filt <- dplyr::bind_rows(references, 
                                         dplyr::select(studies_filt, 
                                                       studyId, ref_complete, DOI))
      # get edges
      edges_filt <- edges %>% dplyr::filter(sampleId %in% samples_filt$sampleId)
      
      # get relative abundance of taxa
      taxon_ab <- edges_filt %>% group_by(taxonId) %>%
        summarise(ab_sum = sum(weight)) %>% 
        ungroup() %>%
        mutate(ab_sum = ab_sum/sum(ab_sum))
      
      # get taxa from edges
      taxonId_filt <- edges_filt %>% dplyr::distinct(taxonId) %>% arrange(.) %>% pull(.)
      
      # filter taxa
      taxa_filt <- taxa %>% dplyr::filter(taxonId %in% taxonId_filt) %>%
        left_join(.,taxon_ab)
      
      # building a list with essential elements ---------------------------------
      # initially sample aggregation is sample and taxonomic aggregation is species
      
      FMBN_list <- list(studies = studies_filt,
                        samples = samples_filt,
                        edges = edges_filt,
                        taxa = taxa_filt,
                        references = reference_filt,
                        version = paste0("This dataset includes ", 
                                         nrow(samples_filt), " samples from FoodMicrobionet, extracted on ",
                                         today(),
                                         ". ",
                                         version),
                        sample_Ids = pull(f_fsamples()[input$ssamplestable_rows_all,2]),
                        sample_agg = "sample",
                        tax_agg = "species"
      )
      
      return(FMBN_list)
    }
  )
  
  # building the reactive aggregate list ------------------------------------
  # will be used in any case to export data

  FMBNexp <- reactive(
    {
      taxa_redux <- switch(input$tax_agg_level,
                           "species" = {FMBNfilt()[[4]] %>% 
                               mutate(a_label = label, t_label = genus)},
                           "genus" = {FMBNfilt()[[4]] %>% 
                               dplyr::select(taxonId:taxonomy, idelevel) %>%
                               mutate(species = NA)
                           },
                           "family" = {FMBNfilt()[[4]] %>% 
                               dplyr::select(taxonId:taxonomy, idelevel) %>%
                               mutate(genus = NA, species = NA)
                           },
                           "class" = {FMBNfilt()[[4]] %>% 
                               dplyr::select(taxonId:taxonomy, idelevel) %>%
                               mutate(order = NA, family = NA, genus = NA, species = NA)
                           }
      )
      if(input$tax_agg_level != "species") {
        taxa_redux <- switch(input$tax_agg_level,
                             "genus" = {taxa_redux %>% 
                                 mutate(id_L6 = str_sub(id_L6, 1, str_locate(id_L6, ";s_")[,1]+3),
                                        taxonomy = str_sub(taxonomy, 1, str_locate(taxonomy, "; s_")[,1]+4),
                                        t_label = genus,
                                        a_label = ifelse(is.na(t_label), label, t_label)
                                 )
                             },
                             "family" = {taxa_redux %>%
                                 mutate(id_L6 = paste0(str_sub(id_L6, 1, str_locate(id_L6, ";g_")[,1]+3),";s__"),
                                        taxonomy = paste0(str_sub(taxonomy, 1, str_locate(taxonomy, "; g_")[,1]+4),"; s__"),
                                        t_label = family,
                                        a_label = ifelse(is.na(t_label), label, t_label)
                                 )
                             },
                             "class" = {taxa_redux %>%
                                 mutate(id_L6 = paste0(str_sub(id_L6, 1, str_locate(id_L6, ";f_")[,1]+3),";g__;s__"),
                                        taxonomy = paste0(str_sub(taxonomy, 1, str_locate(taxonomy, "; f_")[,1]+4),"; g__; s__"),
                                        t_label = class,
                                        a_label = ifelse(is.na(t_label), label, t_label)
                                 )
                             }
        )
      }
      # adding info to the edges
      edges_exp <- left_join(FMBNfilt()[[3]], 
                             dplyr::select(taxa_redux, -id_L6, -idelevel, -label, -taxonomy, -t_label))
      edges_exp <- left_join(edges_exp, 
                             dplyr::select(FMBNfilt()[[2]], sampleId, label_2, llabel, foodId))              
      
      # cast 
      col_agg_var <- switch(input$sample_agg_level,
                            "sample" = "label_2",
                            "exp. code" = "llabel")
      OTUtableagg <- dcast(edges_exp, as.formula(paste("a_label", col_agg_var, sep = "~")), 
                           sum, value.var = "weight", drop = FALSE, fill = 0)
      
      # transforming in proportions
      sumcheck <- colSums(OTUtableagg[2:ncol(OTUtableagg)])
      tOTUm <- t(OTUtableagg[,2:ncol(OTUtableagg)])/sumcheck
      colnames(tOTUm) <- OTUtableagg$a_label
      OTUtableagg[,2:ncol(OTUtableagg)] <- t(tOTUm)
      
      # rebuilding a list of taxa, in case of aggregation
  
      taxa_list <- taxa_redux %>% distinct(a_label, .keep_all = T) %>%
        dplyr::select(-taxonId, -label, -idelevel, -t_label)
      
      
      # getting nreads
      nreads <- switch(col_agg_var,
                       "label_2" = FMBNfilt()[[2]] %>% 
                         dplyr::select(label = label_2, n_reads2) %>% arrange(label),
                       "llabel" = FMBNfilt()[[2]] %>%
                         dplyr::select(label = llabel, n_reads2) %>% arrange(label) %>%
                         mutate(n_reads2 = median(n_reads2)) %>% 
                         group_by(label) %>%
                         summarise(n_reads2 = median(n_reads2)*n())
      )
      # transforming tOTUm
      tOTUm <- round(tOTUm * nreads$n_reads2)
      # prepping the sample table
      sample_table <- switch(input$sample_agg_level,
                             "sample" = {
                               FMBNfilt()[[2]] %>% 
                                 dplyr::select(studyId, label = label_2, llabel, s_type, 
                                               n_reads2:SRA_run)
                               },
                             "exp. code" = {
                               FMBNfilt()[[2]] %>%
                                 dplyr::select(-n_reads2) %>%
                                 mutate(studyId = NA, label = llabel, description = " ", 
                                        target1 = NA, target2 = NA) %>%
                                 distinct(label, .keep_all = TRUE) %>%
                                 left_join(., nreads, by = "label") %>%
                                 dplyr::select(studyId, label, llabel, s_type, n_reads2, 
                                        foodId:SRA_run)
                               }
                             )
      sample_table <- as.data.frame(sample_table)
      rownames(sample_table) <- sample_table$label
      
      # making an object for import in Gephi  or Cytoscape ---------------------------
      
      edges <- OTUtableagg %>% tidyr::gather(key = "Source", value = "weight", -a_label) %>%
        dplyr::mutate(weight = weight*100) %>% 
        dplyr::select(Source, Target = a_label, weight) %>%
        dplyr::filter(weight > 0)
      
      sample_nodes <- switch(input$sample_agg_level,
                             "sample" = {sample_table %>% 
                                 dplyr::select(label, node_type = s_type, L1:L6, 
                                               info1 = llabel, info2 = foodId, info3 = nature,
                                               info4 = process, info5 = spoilage, 
                                               info6 = description, info7 = target1,
                                               info8 = target2)
                             },
                             "exp. code" = {sample_table %>%
                                 mutate(info1 = label, 
                                        info6 = as.character(NA), 
                                        info7 = as.character(NA), 
                                        info8 = as.character(NA)) %>%
                                 dplyr::select(label, node_type = s_type, L1:L6, 
                                               info1, info2 = foodId, info3 = nature,
                                               info4 = process, info5 = spoilage, info6:info8)
                               
                             }
      )
      OTU_nodes <- taxa_list %>% 
        mutate(node_type = "OTU", L1 = phylum,
               L4 = class, L6 = family, info1 = genus,
               info2 = species, info3 = NA, info4 = NA,
               info5 = NA, info6 = id_L6, info7 = NA,
               info8 = NA) %>%
        dplyr::select(label = a_label, node_type, L1:info8)
      
      nodes <- bind_rows(sample_nodes, OTU_nodes)  
      
      # make a igraph object
      FMBNigraph <- graph_from_data_frame(edges, directed = F, vertices = nodes)
      
      # adding abundance to taxa_list
      taxa_ab <- edges %>% group_by(Target) %>%
        summarise(ab_sum = sum(weight)) %>% 
        ungroup() %>%
        mutate(ab_sum = ab_sum/sum(ab_sum)) %>%
        dplyr::rename(a_label = Target)
      taxa_list <- left_join(taxa_list, taxa_ab)
      
      # assemble the list
      
      FMBN_agg <- list(studies = FMBNfilt()[[1]],
                       OTU_table = tOTUm,
                       OTU_table_relf = OTUtableagg,
                       sample_metadata = sample_table,
                       taxa_metadata = taxa_list %>% dplyr::select(label = a_label, 
                                                                   domain:taxonomy, ab_sum),
                       edge_table = edges,
                       node_table = nodes,
                       i_graph = FMBNigraph,
                       references = FMBNfilt()[[5]],
                       version = FMBNfilt()[[6]],
                       sample_agg = input$sample_agg_level,
                       tax_agg = input$tax_agg_level)
      return(FMBN_agg)
      }
  )
  
  
  # the outputs, aggregate tab ----------------------------------------------
  
  # data on selection, aggregate tab
  output$subset_info <- renderText(FMBNfilt()[[6]])
  
  # study table, aggregate tab, not affected by aggregation
  output$f_studies_table <- DT::renderDataTable({
    FMBNfilt()[[1]] %>% 
      mutate(
        target = paste(target, region),
        SRA_ENA = ifelse(is.na(Seq_accn),
                         "not available", 
                         paste0("<a href='",
                                Seq_accn_link,
                                "' target='_blank'>",
                                Seq_accn,"</a>")),
        DOIlink = ifelse(DOI == "unpublished data",
                         "unpublished data",
                         paste0("<a href='",
                                DOI_link,
                                "' target='_blank'>",
                                DOI,"</a>"))
      ) %>%
      dplyr::select(studyId, FMBN_version, target, platform, samples, SRA_ENA,
                    ref_short, DOIlink, short_descr)
    },
    rownames = F,
    escape = c(7,9),
    colnames = c("study" = "studyId", 
                 "version" = "FMBN_version", 
                 "SRA/ENA" = "SRA_ENA",
                 "reference" = "ref_short",
                 "Seq. accn. link" = "SRA_ENA",
                 "DOI link" = "DOIlink",
                 "description" = "short_descr"),
    options = list(pageLength = 3,
                   lengthMenu = c(3, 5, 10),
                   columnDefs = list(
                     list(
                       targets = 8,
                       render = JS(
                         "function(data, type, row, meta) {",
                         "return type === 'display' && data.length > 50 ?",
                         "'<span title=\"' + data + '\">' + data.substr(0, 50) + '...</span>' : data;","}"
                         )
                       )
                     )
                   )
    )
  # sample table, aggregate tab, affected by sample_agg_level and show_what
  
  output$f_samples_table <- DT::renderDataTable(
    if(input$sample_agg_level == "samples" | input$show_what == "original"){
      dplyr::select(FMBNfilt()[[2]],
             studyId, label = label_2, s_type, n_reads2, foodId, llabel, L1, L6,
             description, nature:target2)
      } else {
        dplyr::select(FMBNexp()[[4]],
               studyId, label, s_type, n_reads2, foodId, llabel, L1, L6,
               description, nature:target2)  
      },
    rownames = F,
    escape = T,
    colnames = c("study" = "studyId",
                 "label" = "label",
                 "sample type" = "s_type", 
                 "reads" = "n_reads2",
                 "code" = "foodId", 
                 "exp. code" = "llabel",
                 "L1 food group" = "L1",
                 "L6 food subgroup" = "L6",
                 "description" = "description",
                 "DNA/cDNA" = "target1",
                 "region" = "target2"),
    options = list(pageLength = 3,
                   lengthMenu = c(3, 5, 10, 25, 20),
                   columnDefs = list(
                     # limits display length of columns 7-9
                     list(
                       targets = c(7,8,9),
                       render = JS(
                         "function(data, type, row, meta) {",
                         "return type === 'display' && data.length > 20 ?",
                         "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;","}"
                       )
                     )
                   )
                   
    )
  )

  # taxa table, aggregate tab, affected by sample_agg_level and show_what

  output$f_taxa_table <- DT::renderDataTable(
    if(input$tax_agg_level == "species" | input$show_what == "original"){
      FMBNfilt()[[4]] %>%
        dplyr::arrange(-ab_sum) %>%
        dplyr::select(label, ab_sum, domain:species,
                      NCBI_outlink, BacterioNet_outlink) %>%
        dplyr::mutate(
          ab_sum = round(ab_sum, digits = 4),
          NCBI_outlink = paste0("<a href='",
                                NCBI_outlink,
                                "' target='_blank'>",
                                NCBI_outlink,"</a>"),
          BacterioNet_outlink = paste0("<a href='",
                                       BacterioNet_outlink,
                                       "' target='_blank'>",
                                       BacterioNet_outlink,"</a>")
        )
    } else {
      FMBNexp()[[5]] %>%
        dplyr::arrange(-ab_sum) %>%
        dplyr::select(label, ab_sum, domain:species) %>%
        dplyr::mutate(
          ab_sum = round(ab_sum, digits = 4),
          NCBI_outlink = ifelse(label == "Other", 
                                "<a href='http://www.foodmicrobionet.org/?p=235' target='_blank'>http://www.foodmicrobionet.org/?p=235</a>",
                                paste0("<a href='",
                                       "http://www.ncbi.nlm.nih.gov/taxonomy/?term=",
                                       label,
                                       "' target='_blank'>",
                                       "http://www.ncbi.nlm.nih.gov/taxonomy/?term=",
                                       label,
                                       "</a>")),
          BacterioNet_outlink = ifelse(label == "Other", 
                                       "<a href='http://www.foodmicrobionet.org/?p=235' target='_blank'>http://www.foodmicrobionet.org/?p=235</a>",
                                       paste0("<a href='",
                                              "http://www.bacterio.net/",
                                              str_to_lower(label),
                                              ".html",
                                              "' target='_blank'>",
                                              "http://www.bacterio.net/",
                                              str_to_lower(label),
                                              ".html",
                                              "</a>"))
          )
    },
    rownames = F,
    escape = c(7,8),
    colnames = c("label" = "label",
                 "rel. abund." = "ab_sum",
                 "domain" = "domain",
                 "phylum" = "phylum",
                 "class" = "class",
                 "family" = "family",
                 "genus" = "genus",
                 "species" = "species",
                 "NCBI" = "NCBI_outlink",
                 "bacterio.net" = "BacterioNet_outlink"),
    options = list(pageLength = 3,
                   lengthMenu = c(3, 5, 10, 25, 20)
    )
  )
  # reference list, aggregate tab
  output$f_refs <- DT::renderDataTable(
    {FMBNfilt()[[5]] %>% mutate(
      DOIlink = ifelse(DOI == "unpublished data",
                       "unpublished data",
                       paste0("<a href='",
                              paste0("https://doi.org/", DOI),
                              "' target='_blank'>",
                              DOI,"</a>")
                       )
      ) %>%
        dplyr::select(-DOI)
      },
    rownames = F,
    escape = c(2),
    colnames = c("Set" = "studyId",
                 "Reference." = "ref_complete",
                 "DOI" = "DOIlink"
                 ),
    options = list(pageLength = 5,
                   lengthMenu = c(5, 10, 25, 20)
    )
  )
  
  # functions, export tab ---------------------------------------------------
  
  # the summary info
  output$mysummary <- DT::renderDataTable(
    data.frame(
      Summary = c("Studies", "Samples", "Food groups", "Food codes", "exp. food codes", "taxa"),
      number= c(nrow(FMBNfilt()[[1]]),
                nrow(FMBNfilt()[[2]]),
                n_distinct(FMBNfilt()[[2]]$L1),
                n_distinct(FMBNfilt()[[2]]$foodId),
                n_distinct(FMBNfilt()[[2]]$llabel),
                nrow(FMBNfilt()[[4]]))
    ),
    rownames = F,
    options = list(
      dom = "t",
      autoWidth = TRUE,
      columnDefs = list(list(width = '10%', targets = c(1, 1)))
      )
  )
  
  output$agg_info <- renderText(paste("The sample aggregation level is ",
                                        input$sample_agg_level, ". ",
                                        "The aggregation level for taxa is ",
                                        input$tax_agg_level, ". There are ",
                                        nrow(FMBNexp()[[3]]), " taxa and ",
                                        ncol(FMBNexp()[[3]]), " samples in your",
                                        " dataset after aggregation.",
                                        sep = "")
                                 )
  
  output$OTU_table_agg <- DT::renderDataTable(
    {OTUdisp <- FMBNexp()[[3]] %>%
      mutate_if(is.numeric, round, digits = 4)
    },
    rownames = F,
    options = list(
      dom = "lftip",
      scrollX = T
    )
  )
  # save selection no aggregation
  observeEvent(input$miniFMBN,
               {
                 saveRDS(FMBNfilt(),
                         file = file.path("output", "minifmbn", 
                                          paste(input$fn_prefix,input$miniFMBN,".RDS",sep=""))
                 )
                 showNotification("urFMBN file saved to /output/minifmbn folder", 
                                  type = "message", duration = 5, closeButton = TRUE)
               })
  
  # save selection after aggregation
  observeEvent(input$aggFMBN,
               {
                 saveRDS(FMBNexp(),
                         file = file.path("output", "aggdata", 
                                          paste(input$fn_prefix,input$aggFMBN,"_agg.RDS",sep=""))
                 )
                 showNotification("aggFMBN file saved to /output/aggdata folder", 
                                  type = "message", duration = 5, closeButton = TRUE)
               })
  
  # create and save save a phyloseq object after aggregation
  observeEvent(input$phyloseq,
               {
                 taxtable <- FMBNexp()[[5]] %>% dplyr::select(label, domain:species)
                 taxtable <- as.data.frame(taxtable)
                 rownames(taxtable) <- taxtable$label
                 taxtable <- as.matrix(taxtable[,2:ncol(taxtable)])
                 # create and a phyloseq class object
                 physeqdata <- phyloseq(otu_table(t(FMBNexp()[[2]]), taxa_are_rows = T), 
                                        tax_table(taxtable), sample_data(FMBNexp()[[4]]))
                 # phyloseq object
                 save(physeqdata,
                      file = file.path("output", "phyloseq", paste(input$fn_prefix,
                                                                   input$phyloseq,
                                                                   "_physeq.Rdata",sep=""))
                 )
                 showNotification("phyloseq file saved to /output/phyloseq folder", 
                                  type = "message", duration = 5, closeButton = TRUE)
               })
  
  # create and save a OTU and sample table for CoNet
  observeEvent(input$conet,
               {
                 OTU_table_CoNet <- as.data.frame(t(FMBNexp()[[2]])) %>% 
                   tibble::rownames_to_column(var = "OTUID") %>%
                   inner_join(., dplyr::select(FMBNexp()[[5]], taxonomy, label), by = c("OTUID" = "label"))
                 write_tsv(OTU_table_CoNet,
                           file.path("output", "conet", paste(input$fn_prefix,
                                                              input$conet,
                                                              "_conet_OTU.txt",sep=""))
                 )
                 write_tsv(remove_rownames(FMBNexp()[[4]]) %>% 
                             dplyr::rename(SampleID = label) %>%
                             dplyr::select(SampleID, setdiff(names(FMBNexp()[[4]]),"label")),
                           file.path("output", "conet", paste(input$fn_prefix,
                                                              input$conet,
                                                              "_conet_samples.txt",sep=""))
                 )
                 showNotification("OTU and sample tables for CoNet saved to /output/conet folder", 
                                  type = "message", duration = 5, closeButton = TRUE)
               }
               )
  
  # create and save a .gml graph for use in Gephi and Cytoscape
  observeEvent(input$gml,
               {
                 write.graph(FMBNexp()[[8]],
                             file.path("output", "gml", paste(input$fn_prefix,
                                                              input$gml,
                                                              "_igraph.gml",sep="")), 
                             format = "gml")
                 showNotification(".gml network saved in /output/gml folder", 
                                  type = "message", duration = 5, closeButton = TRUE)
               }
  )
  
}


## questo per caricare eventuali helpers, che devono essere nella stessa cartella
# source("helpers.R")

# Create Shiny app ----
shinyApp(ui = ui, server = server)
