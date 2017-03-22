library(shinydashboard)
library(shiny)
library(d3heatmap)
library(pcaExplorer)
library(DESeq2)
library(ggplot2)
library(shinyAce)
library(DT)
library(knitr)
library(rmarkdown)
library(pheatmap)
library(topGO)
library(UpSetR)
library(BiocParallel)
library(rintrojs)
library(IHW)
library(goseq)



modes <- shinyAce::getAceModes()
themes <- shinyAce::getAceThemes()

ideal_env <<- new.env(parent = emptyenv())


## upload max 300mb files - can be changed if necessary
options(shiny.maxRequestSize=300*1024^2)


footer <- function(){
  tags$div(
    class = "panel-footer",
    style = "text-align:center",
    tags$div(
      class = "foot-inner",
      list(
        # hr(),
        "ideal is a project developed by Federico Marini in the Bioinformatics division of the ",
        tags$a(href="http://www.unimedizin-mainz.de/imbei","IMBEI"),
        "- Institute for Medical Biostatistics, Epidemiology and Informatics",br(),
        "License: ",tags$a(href="https://opensource.org/licenses/MIT","MIT"), br(),

        "Development of the ideal package is on ",
        tags$a(href="https://github.com/federicomarini/ideal", "GitHub")
      )
    )
  )
}


ideal_ui <- shinydashboard::dashboardPage(
  shinydashboard::dashboardHeader(
    title = paste0("ideal - Interactive Differential Expression AnaLysis ",
                   packageVersion("ideal")),
    titleWidth = 900,

    # TODO:
    # http://stackoverflow.com/questions/31440564/adding-a-company-logo-to-shinydashboard-header
    # replace text with image
    # ideal_header$children[[2]]$children <- tags$a(href='https://github.com/federicomarini/ideal',
    # tags$img(src='ideal_logo_v2.png',height='50',width='200'))
    # title = tags$a(href='https://github.com/federicomarini/ideal',
    #                tags$img(src='ideal_logo_v2.png',height='50',width='200')),

    # task menu for saving state to environment or binary data
    shinydashboard::dropdownMenu(type = "tasks",icon = icon("cog"),badgeStatus = NULL, # something to change the top message? maybe file an issue @shinydashboard development
                                 notificationItem(
                                   text = actionButton("task_exit_and_save","Exit ideal & save",
                                                       class = "btn_no_border",
                                                       onclick = "setTimeout(function(){window.close();}, 100); "),
                                   icon = icon("sign-out"),status = "primary"),
                                 menuItem(
                                   text = downloadButton("task_state_save","Save State as .RData"))
    )
  ),


  dashboardSidebar(
    width = 280,
    menuItem("App settings",icon = icon("cogs"),
             uiOutput("color_by"),
             uiOutput("available_genes"),

             #
             #              selectizeInput(
             #                inputId = 'available_genes', label = 'Select Something',
             #                choices = NULL,
             #                multiple = TRUE,
             #                selected = 1
             #              ),

             numericInput("FDR","False Discovery Rate",value = 0.05, min = 0, max = 1, step = 0.01)

    ),
    menuItem("Plot export settings", icon = icon("paint-brush")),
    menuItem("Quick viewer", icon = icon("flash"),
             fluidRow(
               fluidRow(column(6,p("Count matrix")), column(6,uiOutput("ok_cm"))),
               fluidRow(column(6,p("Experimental design")), column(6,uiOutput("ok_ed"))),
               fluidRow(column(6,p("DESeqDataset")), column(6,uiOutput("ok_dds"))),
               fluidRow(column(6,p("Annotation")), column(6,uiOutput("ok_anno"))),
               fluidRow(column(6,p("Results")), column(6,uiOutput("ok_resu")))
             )),
    menuItem("First steps help", icon = icon("question-circle"),
             actionButton("btn", "Click me for a quick tour", icon("info"),
                          style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4")
    )
  ),



  dashboardBody(
    introjsUI(),
    # must include in UI



    ## Define output size and style of error messages, and also the style of the icons e.g. check
    ## plus, define the myscrollbox div to prevent y overflow when page fills up
    tags$head(
      tags$style(HTML("
                      .shiny-output-error-validation {
                      font-size: 15px;
                      color: forestgreen;
                      text-align: center;
                      }
                      .icon-done {
                      color: green;
                      }
                      #myScrollBox{
                      overflow-y: scroll;

                      .dataTables_wrapper{
                      overflow-x: scroll;
                      }

                      }
                      "))
      ),

    # value boxes to always have an overview on the available data
    fluidRow(
      valueBoxOutput("box_ddsobj"),
      valueBoxOutput("box_annobj"),
      valueBoxOutput("box_resobj")
    ),

    ## main structure of the body for the dashboard
    div(
      id = "myScrollBox", # trick to have the y direction scrollable
      tabBox(
        width=12,


        tabPanel(
          # "Welcome!",  icon = icon("info-circle"),
          title = "Welcome!",  icon = icon("home"), value="tab-welcome",
          # carouselPanel(
          #   img(src = "www/ideal_logo_v2.png"),
          #   img(src = "ideal_logo_v2.png"),
          #   img(src = "ideal_logo_v2.png")
          # ),


          ### TODO: proof of principle it works with the carousel, to display functionality at once
          # bs_carousel(id = "the_beatles", use_indicators = TRUE) %>%
          #   bs_append(
          #     content = bs_carousel_image(src = image_uri("inst/extdata/ideal_logo_v2.png")),
          #     caption = bs_carousel_caption("John Lennon", "Rhythm guitar, vocals")
          #   ) %>%
          #   bs_append(
          #     content = bs_carousel_image(src = image_uri("figure/unnamed-chunk-4-1.png")),
          #     caption = bs_carousel_caption("Paul McCartney", "Bass guitar, vocals")
          #   ),



          ## TODO: explore possibility to put a carousel of images: https://github.com/dcurrier/carouselPanel/
          # carouselPanel(
          #   plotOutput("distPlot1"),
          #   plotOutput("distPlot2")
          # ),
          # img(src = "ideal_logo_v2.png"),


          fluidRow(
            column(
              width = 8,
              introBox(includeMarkdown(system.file("extdata", "welcome.md",package = "ideal"))
                       ,data.step = 1,data.intro = "welcome on board!"),

              br(),br(),

              p("If you see a grey box like this one open below..."),

              shinyBS::bsCollapse(id = "help_welcome",open = "Help", # alt: "Help"
                                  # think of a general trigger for this? something like a variable that assumes NULL
                                  shinyBS::bsCollapsePanel("Help", includeMarkdown(system.file("extdata", "help_welcome.md",package = "ideal")))
              ),

              actionButton("introexample", "If you see a button like this...", icon("info"),
                           style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"),
              p("... you can click on that to start a tour based on introJS"),
              br(),br(),

              introBox(includeMarkdown(system.file("extdata", "instructions.md",package = "ideal")),data.step = 2,data.intro = "follow here")
            )
          )

        ),


        tabPanel(
          "Data Setup",icon = icon("upload"), # value="tab-ds",
          value = "tab-datasetup",

          headerPanel("Setup your data for the analysis"),

          fluidRow(
            column(
              width = 8,
              shinyBS::bsCollapse(id = "help_datasetup",open = NULL, # alt: "Help"
                                  # think of a general trigger for this? something like a variable that assumes NULL
                                  shinyBS::bsCollapsePanel("Help",includeMarkdown(system.file("extdata", "help_datasetup.md",package = "ideal")))
              )
            )
          ),

          actionButton("tour_datasetup", "Click me for a quick tour of the section", icon("info"),
                       style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"), br(),

          introBox(box(width = 12, title = "Step 1", status = "danger", solidHeader = TRUE,
                       h2("Upload your count matrix and the info on the experimental design"),

                       fluidRow(
                         column(
                           width = 4,
                           uiOutput("upload_count_matrix"),
                           uiOutput("upload_metadata"),
                           br(),
                           "... or you can also ",
                           actionButton("btn_loaddemo", "Load the demo airway data", icon = icon("play-circle"),
                                        class = "btn btn-info"),br(), p()
                         )
                       ),

                       fluidRow(
                         column(
                           width = 6,
                           box(width = NULL, title = "Count matrix preview",status = "primary",
                               solidHeader = TRUE,collapsible = TRUE, collapsed = TRUE,
                               fluidRow(
                                 column(
                                   width = 12,
                                   offset = 0.5,
                                   DT::dataTableOutput("dt_cm"))
                               )
                           )
                         ),
                         column(
                           width = 6,
                           box(width = NULL, title = "Experimental design preview",status = "primary",
                               solidHeader = TRUE,collapsible = TRUE, collapsed = TRUE,
                               fluidRow(
                                 column(
                                   width = 12,
                                   offset = 0.5,
                                   DT::dataTableOutput("dt_ed"))
                               )
                           )
                         )
                       )


          ), data.step = 3,data.intro = "upload your data and do stuff"),
          # h2("Step 1: Upload your count matrix and the info on the experimental design"),


          uiOutput("ui_step2"),

          # hr(),
          fluidRow(
            column(
              width = 6,
              uiOutput("ui_stepanno")
              ## this ideally populates also the list of genes of interest to choose among
            ),
            column(
              width = 6,
              uiOutput("ui_stepoutlier")
            )
          ),
          # hr(),

          uiOutput("ui_step3")
        ),




        tabPanel(
          "Counts Overview",
          icon = icon("eye"),
          conditionalPanel(
            condition="!output.checkdds",

            headerPanel("Get an overview on your data"),

            fluidRow(
              column(
                width = 8,
                shinyBS::bsCollapse(id = "help_countsoverview",open = NULL, # alt: "Help"
                                    # think of a general trigger for this? something like a variable that assumes NULL
                                    shinyBS::bsCollapsePanel("Help",includeMarkdown(system.file("extdata", "help_overview.md",package = "ideal")))
                )
              )
            ),

            actionButton("tour_countsoverview", "Click me for a quick tour of the section", icon("info"),
                         style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"), br(),


            ### to control colors of action buttons
            # actionButton("run", "Run Analysis", icon("paper-plane"),
            # style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),

            selectInput("countstable_unit", label = "Data scale in the table",
                        choices = list("Counts (raw)" = "raw_counts",
                                       "Counts (normalized)" = "normalized_counts",
                                       "Regularized logarithm transformed" = "rlog_counts",
                                       "Log10 (pseudocount of 1 added)" = "log10_counts",
                                       "TPM (Transcripts Per Million)" = "tpm_counts")),

            DT::dataTableOutput("showcountmat"),
            downloadButton("downloadData","Download", class = "btn btn-success"),
            hr(),
            fluidRow(
              column(
                width = 8,
                h3("Basic summary for the counts"),
                p("Number of uniquely aligned reads assigned to each sample"),
                # verbatimTextOutput("reads_summary"),
                wellPanel(
                  fluidRow(
                    column(
                      width = 6,
                      numericInput("threshold_rowsums","Threshold on the row sums of the counts",value = 0, min = 0)),
                    column(
                      width = 6,
                      numericInput("threshold_rowmeans","Threshold on the row means of the normalized counts",value = 0, min = 0))
                  )),
                p("According to the selected filtering criteria, this is an overview on the provided count data"),
                verbatimTextOutput("detected_genes"),

                selectInput("filter_crit",label = "Choose the filtering criterium",
                            choices = c("row means", "row sums"), selected = "row means"),

                actionButton("featfilt_dds", "Filter the DDS object",class = "btn btn-primary")
              )
            ),


            h3("Sample to sample scatter plots"),
            selectInput("corr_method","Correlation method palette",choices = list("pearson","spearman")),
            p("Compute sample to sample correlations on the normalized counts - warning, it can take a while to plot all points (depending mostly on the number of samples you provided)."),
            actionButton("compute_pairwisecorr", "Run", class = "btn btn-primary"),
            uiOutput("pairwise_plotUI"),
            uiOutput("heatcorr_plotUI")

          ),
          conditionalPanel(
            condition="output.checkdds",
            h2("You did not create the dds object yet. Please go the main tab and generate it")
          )
        ),



        tabPanel(
          "Extract Results", icon = icon("table"),

          # see: http://stackoverflow.com/questions/21609436/r-shiny-conditionalpanel-output-value?noredirect=1&lq=1
          conditionalPanel(
            condition="!output.checkdds",

            headerPanel("Extract and inspect the DE results"),

            fluidRow(
              column(
                width = 8,
                shinyBS::bsCollapse(id = "help_extractresults",open = NULL, # alt: "Help"
                                    # think of a general trigger for this? something like a variable that assumes NULL
                                    shinyBS::bsCollapsePanel("Help",includeMarkdown(system.file("extdata", "help_results.md",package = "ideal")))
                )
              )
            ),

            actionButton("tour_results", "Click me for a quick tour of the section", icon("info"),
                         style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"), br(),



            # conditionalPanel(
            #   condition="output.checkresu==0",
            #   h2('RESU not provided')
            # ),

            # "Data Overview", icon = icon("eye"),

            fluidRow(
              column(
                width = 6,
                uiOutput("choose_fac")
              )
            ),
            fluidRow(
              column(
                width = 4,
                # factor as covariate
                wellPanel(
                  width = 4,
                  uiOutput("fac1"),
                  uiOutput("fac2"),
                  # continuous covariate
                  uiOutput("facnum")
                )

              ),
              column(
                width = 4,
                # factor with > 2 levels
                wellPanel(
                  width = 4,
                  uiOutput("lrtavailable"),
                  uiOutput("lrtfull"),
                  uiOutput("lrtreduced")
                ),

                uiOutput("runlrt")
              )

            ),

            ## general options for result function
            # alpha is set via FDR on the left side
            fluidRow(
              column(
                width = 4,
                wellPanel(
                  selectInput("resu_indfil",label = "Apply independent filtering automatically",
                              choices = c(TRUE,FALSE), selected = TRUE),
                  selectInput("resu_addmle",label = "Add the unshrunken MLE of log2 fold change",
                              choices = c(TRUE,FALSE), selected = TRUE),
                  selectInput("resu_ihw", "Use Independent Hypothesis Weighting (IHW) as a filtering function",
                              choices = c(TRUE, FALSE), selected = FALSE)
                )
              )
            ),
            #, evtl also the *filter* parameter of the function, i.e. baseMean if not specified

            fluidRow(
              column(
                width = 6,
                uiOutput("runresults"),
                uiOutput("store_result"),
                verbatimTextOutput("diyres_summary")
              )
            ),





            DT::dataTableOutput("table_res"),
            fluidRow(
              column(
                width = 6,
                plotOutput("pvals_hist")
              ),
              column(
                width = 6,
                plotOutput("logfc_hist")
              )
            )
          ),
          conditionalPanel(
            condition="output.checkdds",
            h2("You did not create the dds object yet. Please go the main tab and generate it")
          )

        ),


        tabPanel(
          "Summary Plots", icon = icon("photo"),
          conditionalPanel(
            condition="!output.checkresu",

            headerPanel("Interactive graphical exploration of the results"),

            fluidRow(
              column(
                width = 8,
                shinyBS::bsCollapse(id = "help_summaryplots",open = NULL, # alt: "Help"
                                    # think of a general trigger for this? something like a variable that assumes NULL
                                    shinyBS::bsCollapsePanel("Help",includeMarkdown(system.file("extdata", "help_plots.md",package = "ideal")))
                )
              )
            ),

            actionButton("tour_plots", "Click me for a quick tour of the section", icon("info"),
                         style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"), br(),

            fluidRow(column(6,
                            h4("MA plot - Interactive!"),
                            plotOutput('plotma', brush = 'ma_brush')),
                     column(6,
                            h4("Zoomed section"),
                            plotOutput("mazoom",click= 'mazoom_click'))
                     # ,
                     # column(4,
                     #        h4("Boxplot for the selected gene"),
                     #        plotOutput("geneplot")
                     # )
            ),
            fluidRow(column(6,
                            h4("Selected gene"),
                            checkboxInput("ylimZero_genes","Set y axis limit to 0",value=TRUE),
                            plotOutput("genefinder_plot")
            ),
            column(6,
                   h4("Gene infobox"),
                   htmlOutput("rentrez_infobox"))

            ),


            fluidRow(column(6,
                            h4("volcano plot"),
                            plotOutput("volcanoplot")
            )),


            # plotOutput("volcanoplot"),
            fluidRow(radioButtons("heatmap_colv","Cluster samples",choices = list("Yes"=TRUE,"No"=FALSE),selected = TRUE)),
            fluidRow(
              column(4,
                     checkboxInput("rowscale",label = "Scale by rows",value = TRUE)),
              column(4,
                     checkboxInput("pseudocounts","use log2(1+counts)",value = TRUE))
            ),
            fluidRow(
              column(6,
                     plotOutput("heatbrush")
              ),
              column(6,
                     d3heatmapOutput("heatbrushD3"))
            ),


            box(
              title = "Brushed table", status = "primary", solidHeader = TRUE,
              collapsible = TRUE, collapsed = TRUE, width = 12,
              fluidRow(DT::dataTableOutput("ma_brush_out")))
          ),
          conditionalPanel(
            condition="output.checkresu",
            h2("You did not create the result object yet. Please go the dedicated tab and generate it")
          )
        ),
        tabPanel(
          "Gene Finder", icon = icon("crosshairs"),
          conditionalPanel(
            condition="!output.checkdds",

            headerPanel("Find your gene(s) of interest"),

            fluidRow(
              column(
                width = 8,
                shinyBS::bsCollapse(id = "help_genefinder",open = NULL, # alt: "Help"
                                    # think of a general trigger for this? something like a variable that assumes NULL
                                    shinyBS::bsCollapsePanel("Help",includeMarkdown(system.file("extdata", "help_genefinder.md",package = "ideal")))
                )
              )
            ),
            actionButton("tour_genefinder", "Click me for a quick tour of the section", icon("info"),
                         style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"), br(),

            fluidRow(
              column(6,checkboxInput("ylimZero_genefinder","Set y axis limit to 0",value=TRUE))),
            fluidRow(
              column(6,
                     plotOutput("bp1")
              ),
              column(6,
                     plotOutput("bp2"))
            ),
            fluidRow(
              column(6,
                     plotOutput("bp3")
              ),
              column(6,
                     plotOutput("bp4"))
            ),

            plotOutput("ma_highlight"),
            DT::dataTableOutput("table_combi"),


            fileInput(inputId = "gl_ma",
                      label = "Upload a gene list file",
                      accept = c("text/csv", "text/comma-separated-values",
                                 "text/tab-separated-values", "text/plain",
                                 ".csv", ".tsv"), multiple = FALSE),
            plotOutput("ma_hl_list"),
            DT::dataTableOutput("table_combi_list")


          ),
          conditionalPanel(
            condition="output.checkdds",
            h2("You did not create the dds object yet. Please go the main tab and generate it")
          )
        ),

        tabPanel(
          "Functional Analysis", icon = icon("list-alt"),
          conditionalPanel(
            condition="!output.checkresu",

            headerPanel("Find functions enriched in gene sets"),

            fluidRow(
              column(
                width = 8,
                shinyBS::bsCollapse(id = "help_functionalanalysis",open = NULL, # alt: "Help"
                                    # think of a general trigger for this? something like a variable that assumes NULL
                                    shinyBS::bsCollapsePanel("Help",includeMarkdown(system.file("extdata", "help_funcanalysis.md",package = "ideal")))
                )
              )
            ),
            actionButton("tour_funcanalysis", "Click me for a quick tour of the section", icon("info"),
                         style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"), br(),

            selectInput("go_cats",label = "Select the GO category(ies) of interest",
                        choices = list("GO Biological Process" = "BP", "GO Molecular Function" = "MF", "GO Cellular Component" = "CC"),
                        selected = "BP",multiple = TRUE
            ),


            tabBox(
              width = NULL,
              id="gse_tabbox",
              tabPanel("UPregu", icon = icon("arrow-circle-up"),
                       fluidRow(column(width = 6,actionButton("button_enrUP", "Perform gene set enrichment analysis on the upregulated genes",class = "btn btn-primary"))),
                       fluidRow(column(width = 6,actionButton("button_enrUP_goseq", "Perform gene set enrichment analysis on the upregulated genes - goseq",class = "btn btn-primary"))),
                       fluidRow(column(width = 6,actionButton("button_enrUP_topgo", "Perform gene set enrichment analysis on the upregulated genes - topGO",class = "btn btn-primary"))),
                       DT::dataTableOutput("DT_gse_up"),
                       DT::dataTableOutput("DT_gse_up_goseq"),
                       fluidRow(
                         column(width = 9, DT::dataTableOutput("DT_gse_up_topgo")),
                         column(width = 3, plotOutput("goterm_heatmap_up_topgo"))
                       )
              ),
              tabPanel("DOWNregu", icon = icon("arrow-circle-down"),
                       fluidRow(column(width = 6,actionButton("button_enrDOWN", "Perform gene set enrichment analysis on the downregulated genes",class = "btn btn-primary"))),
                       fluidRow(column(width = 6,actionButton("button_enrDOWN_goseq", "Perform gene set enrichment analysis on the downregulated genes - goseq",class = "btn btn-primary"))),
                       fluidRow(column(width = 6,actionButton("button_enrDOWN_topgo", "Perform gene set enrichment analysis on the downregulated genes - topGO",class = "btn btn-primary"))),
                       DT::dataTableOutput("DT_gse_down"),
                       DT::dataTableOutput("DT_gse_down_goseq"),
                       fluidRow(
                         column(width = 9, DT::dataTableOutput("DT_gse_down_topgo")),
                         column(width = 3, plotOutput("goterm_heatmap_down_topgo"))
                       )
              ),
              tabPanel("UPDOWN", icon = icon("arrows-v"),
                       fluidRow(column(width = 6,actionButton("button_enrUPDOWN", "Perform gene set enrichment analysis on the up- and downregulated genes",class = "btn btn-primary"))),
                       fluidRow(column(width = 6,actionButton("button_enrUPDOWN_goseq", "Perform gene set enrichment analysis on the up- and downregulated genes - goseq",class = "btn btn-primary"))),
                       fluidRow(column(width = 6,actionButton("button_enrUPDOWN_topgo", "Perform gene set enrichment analysis on the up- and downregulated genes - topGO",class = "btn btn-primary"))),
                       DT::dataTableOutput("DT_gse_updown"),
                       DT::dataTableOutput("DT_gse_updown_goseq"),
                       fluidRow(
                         column(width = 9, DT::dataTableOutput("DT_gse_updown_topgo")),
                         column(width = 3, plotOutput("goterm_heatmap_updown_topgo"))
                       )
              ),
              tabPanel("List1", icon = icon("list"),
                       fileInput(inputId = "gl1",
                                 label = "Upload a gene list file",
                                 accept = c("text/csv", "text/comma-separated-values",
                                            "text/tab-separated-values", "text/plain",
                                            ".csv", ".tsv"), multiple = FALSE),
                       fluidRow(column(width = 6,actionButton("button_enrLIST1", "Perform gene set enrichment analysis on the genes in list1",class = "btn btn-primary"))),
                       fluidRow(column(width = 6,actionButton("button_enrLIST1_goseq", "Perform gene set enrichment analysis on the list1 genes - goseq",class = "btn btn-primary"))),
                       fluidRow(column(width = 6,actionButton("button_enrLIST1_topgo", "Perform gene set enrichment analysis on the list1 genes - topGO",class = "btn btn-primary"))),
                       DT::dataTableOutput("DT_gse_list1"),
                       DT::dataTableOutput("DT_gse_list1_goseq"),
                       fluidRow(
                         column(width = 9, DT::dataTableOutput("DT_gse_list1_topgo")),
                         column(width = 3, plotOutput("goterm_heatmap_l1_topgo"))
                       )

              ),
              tabPanel("List2", icon = icon("list-alt"),
                       fileInput(inputId = "gl2",
                                 label = "Upload a gene list file",
                                 accept = c("text/csv", "text/comma-separated-values",
                                            "text/tab-separated-values", "text/plain",
                                            ".csv", ".tsv"), multiple = FALSE),
                       fluidRow(column(width = 6,actionButton("button_enrLIST2", "Perform gene set enrichment analysis on the genes in list2",class = "btn btn-primary"))),
                       fluidRow(column(width = 6,actionButton("button_enrLIST2_goseq", "Perform gene set enrichment analysis on the list2 genes - goseq",class = "btn btn-primary"))),
                       fluidRow(column(width = 6,actionButton("button_enrLIST2_topgo", "Perform gene set enrichment analysis on the list2 genes - topGO",class = "btn btn-primary"))),
                       DT::dataTableOutput("DT_gse_list2"),
                       DT::dataTableOutput("DT_gse_list2_goseq"),
                       fluidRow(
                         column(width = 9, DT::dataTableOutput("DT_gse_list2_topgo")),
                         column(width = 3, plotOutput("goterm_heatmap_l2_topgo"))
                       )
              )
            ),







            ## will put collapsible list elements? or multi tab panel? or something to select on the left, and operate output-wise on the right e.g. venn diagrams or table for gene set enrichment
            # h3("custom list 3 - handpicked") # use the select input from the left column?
            # ,verbatimTextOutput("debuggls"),

            # verbatimTextOutput("printUPgenes"),
            # verbatimTextOutput("debuglists"),

            h2("Intersection of gene sets"),

            fluidRow(
              column(width = 4,
                     checkboxInput("toggle_updown","Use up and down regulated genes", TRUE),
                     checkboxInput("toggle_up","Use up regulated genes", FALSE),
                     checkboxInput("toggle_down","Use down regulated genes", FALSE)
              ),
              column(width = 4,
                     checkboxInput("toggle_list1","Use list1 genes", TRUE),
                     checkboxInput("toggle_list2","Use list2 genes", FALSE),
                     checkboxInput("toggle_list3","Use list3 genes", FALSE)
              )
            ),


            fluidRow(
              column(width = 6,plotOutput("vennlists"), offset = 3)),
            fluidRow(
              column(width = 6,plotOutput("upsetLists"), offset = 3))


          ),
          conditionalPanel(
            condition="output.checkresu",
            h2("You did not create the result object yet. Please go the dedicated tab and generate it")
          )
        ),


        tabPanel(
          "Report Editor",
          icon = icon("pencil"),

          headerPanel("Create, view and export a report of your analysis"),

          fluidRow(
            column(
              width = 8,
              shinyBS::bsCollapse(id = "help_reporteditor",open = NULL, # alt: "Help"
                                  # think of a general trigger for this? something like a variable that assumes NULL
                                  shinyBS::bsCollapsePanel("Help",includeMarkdown(system.file("extdata", "help_report.md",package = "ideal")))
              )
            )
          ),

          actionButton("tour_report", "Click me for a quick tour of the section", icon("info"),
                       style="color: #ffffff; background-color: #0092AC; border-color: #2e6da4"), br(),

          fluidRow(
            column(
              width = 6,
              box(
                title = "markdown options", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 9, collapsed = TRUE,
                radioButtons("rmd_dl_format", label = "Choose Format:", c("HTML" = "html", "R Markdown" = "rmd"), inline = TRUE),
                textInput("report_title", "Title: "),
                textInput("report_author", "Author: "),
                radioButtons("report_toc", "Table of Contents", choices = list("Yes" = "true", "No" = "false")),
                radioButtons("report_ns", "Number sections", choices = list("Yes" = "true", "No" = "false")),
                selectInput("report_theme", "Theme", choices = list("Default" = "default", "Cerulean" = "cerulean",
                                                                    "Journal" = "journal", "Flatly" = "flatly",
                                                                    "Readable" = "readable", "Spacelab" = "spacelab",
                                                                    "United" = "united", "Cosmo" = "cosmo")),
                radioButtons("report_echo", "Echo the commands in the output", choices = list("Yes" = "TRUE", "No" = "FALSE")))),
            column(
              width = 6,
              box(
                title = "editor options", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 9, collapsed = TRUE,
                checkboxInput("enableAutocomplete", "Enable AutoComplete", TRUE),
                conditionalPanel(
                  "input.enableAutocomplete",
                  wellPanel(
                    checkboxInput("enableLiveCompletion", "Live auto completion", TRUE),
                    checkboxInput("enableRCompletion", "R code completion", TRUE)
                  )
                ),

                selectInput("mode", "Mode: ", choices=shinyAce::getAceModes(), selected="markdown"),
                selectInput("theme", "Theme: ", choices=shinyAce::getAceThemes(), selected="solarized_light"))
            )
            # ,
            # column( # kept for debugging purposes!
            #   width = 6,
            #   verbatimTextOutput("loadedRmd")
            # )
          ),
          fluidRow(
            column(3,
                   actionButton("updatepreview_button", "Update report",class = "btn btn-primary"),p()
            ),
            column(3, downloadButton("saveRmd", "Generate & Save",class = "btn btn-success"))
          ),

          tabBox(
            width = NULL,
            id="report_tabbox",
            tabPanel("Report preview",
                     icon = icon("file-text"),
                     htmlOutput("knitDoc")
            ),

            tabPanel("Edit report",
                     icon = icon("pencil-square-o"),
                     aceEditor("acereport_rmd", mode="markdown",theme = "solarized_light",autoComplete = "live",
                               value=readLines(system.file("extdata", "irt.Rmd",package = "ideal")),
                               height="800px"))
          )
        ),
        tabPanel(
          "About", icon = icon("institution"),

          # headerPanel("Information on ideal/session"),

          fluidRow(
            column(
              width = 8,
              includeMarkdown(system.file("extdata", "about.md",package = "ideal")),

              verbatimTextOutput("sessioninfo")
            )
          )
        )
        # ,tabPanel(
        #   "devel", icon = icon("github")
        #   # ,
        #   # verbatimTextOutput("debugihw"),
        #   #
        #   # plotOutput("ihwp1"),
        #   # plotOutput("ihwp2"),
        #   # plotOutput("ihwp3"),
        #   # plotOutput("ihwp4"),
        #   # plotOutput("ihwp5"),
        #   # plotOutput("ihwp6"),
        #   # plotOutput("ihwp7"),
        #   # plotOutput("ihwp8")
        #
        # )

      )
    )
    ,footer() # TODO: decide where to place the footer
      ),


  skin="blue"
      )
