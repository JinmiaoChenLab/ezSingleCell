require(shiny)
require(Seurat)
require(dplyr)
require(Matrix)
require(plotly)    ## interactive plots
require(ezSingleCell)
require(DT)    ## customisable datatable output
require(webshot)    ## dependecy for exporting plotly to PDF

shinyUI(fluidPage(
    titlePanel("ezSingleCell: Seurat analysis of scRNAseq data"),
    hr(),

    fluidRow(
        ##------Sidebar---------
        column(3,
               h4('Load Data:'),
               wellPanel(
                   fileInput(inputId = 'tpmFiles',
                             label = "TPM/RPKM file",
                             multiple = FALSE,
                             accept = c("text/csv",
                                        "text/comma-separated-values,text/plain",
                                        ".Robj")),
                   fileInput(inputId = 'cellAnnoFiles',
                             label = "Cell annotation file",
                             multiple = FALSE,
                             accept = c("text/csv",
                                        "text/comma-separated-values,text/plain")),
                   checkboxInput(inputId = 'norm',
                             label = "Normalise?",
                             value = TRUE),
                   fluidRow(
                     column(6,
                            textInput(inputId = "delim",
                                      label = "Name delimiter",
                                      value = "_")
                     ),
                     column(6,
                            numericInput(inputId = "field",
                                         label = "Field",
                                         value = 1,
                                         min = 1)
                     )
                   ),
                   numericInput(inputId = "expThres",
                                label = "Expression cut off",
                                value = 1,
                                min = 0.1),
                   fluidRow(
                       column(6,
                              numericInput(inputId = "min.genes",
                                           label = "Min. genes",
                                           value = 200,
                                           min = 1)
                              ),
                       column(6,
                              numericInput(inputId = "min.cells",
                                           label = "Min. cells",
                                           value = 3,
                                           min = 1)
                              )
                   ),
                   textInput(inputId = "projName",
                             label = "Project Name",
                             value = "Seurat_analysis"),
                   fluidRow(
                     column(6,
                            actionButton("loadButton", "Load Data", icon = icon("hand-o-right"))
                     ),
                     column(6,
                            actionButton("reset", "Reset Data", icon = icon("repeat"))
                     )
                   )
               ),

               ##------Plot download---------
               h4("Export to PDF:"),
               wellPanel(
                 ## Conditional panel for different plots
                 conditionalPanel(" input.QC == 'QC_panel1' && input.tabs == 'QC plots' ",
                                  actionButton("PDFa", "Download Violin Plot in PDF", icon = icon("download"))
                 ),
                 conditionalPanel(" input.QC == 'QC_panel2' && input.tabs == 'QC plots' ",
                                  actionButton("PDFb", "Download Cell Plot in PDF", icon = icon("download"))
                 ),
                 conditionalPanel(" input.tabs == 'Variable Gene Plot' ",
                                  actionButton("PDFc", "Download Variable Genes Plot in PDF", icon = icon("download"))
                 ),
                 conditionalPanel(" input.Pca == 'P_panel1' && input.tabs == 'PCA' ",
                                  actionButton("PDFd", "Download PCA Plot(PDF) and co-ordinates", icon = icon("download"))
                 ),
                 conditionalPanel(" input.Pca == 'P_panel2' && input.tabs == 'PCA' ",
                                  actionButton("PDFe", "Download PC Gene plots and PC table", icon = icon("download"))
                 ),
                 conditionalPanel(" input.diag == 'D_panel1' && input.tabs == 'Diagnostic' ",
                                  actionButton("PDFg", "Download Jackstraw Plot in PDF", icon = icon("download"))
                 ),
                 conditionalPanel(" input.diag == 'D_panel2' && input.tabs == 'Diagnostic' ",
                                  actionButton("PDFh", "Download Elbow Plot in PDF", icon = icon("download"))
                 ),
                 conditionalPanel(" input.tabs == 'TSNE' ",
                                  actionButton("PDFi", "Download TSNE Plot(PDF) and co-ordinates", icon = icon("download"))
                 ),
                 conditionalPanel(" input.tabs == 'DEGs' ",
                                  actionButton("PDFj", "Download DEG results", icon = icon("download"))
                 ),
                 ## ensure no spill over in button text
                 tags$head(
                   tags$style(HTML('
                                   .btn {
                                   white-space: normal;
                                   }'
                                )
                   )
                 ),
                 ## Conditional is separate from pdf options
                 hr(),
                 fluidRow(
                   column(6,
                          sliderInput(inputId="pdf_w", label = "PDF width(in):",
                                      min=3, max=20, value=8, width=100, ticks=F)
                   ),
                   column(6,
                          sliderInput(inputId="pdf_h", label = "PDF height(in):",
                                      min=3, max=20, value=8, width=100, ticks=F)
                   )),

                 actionButton("OpenDir", "Open download folder", icon = icon("folder"))
               ),

               ##------Save Data---------
               hr(),
               actionButton("saveButton", "Save Data", icon = icon("hand-o-right")),

               hr(),
               h4(tags$a(href="mailto:a0124008@u.nus.edu,
                         Chen_Jinmiao@immunol.a-star.edu.sg?subject=[ezSingleCell-question]",
                         "Contact Us")),
               imageOutput("logo", height = "60px")
        ),
        ##------Main area---------
        column(9,
               tabsetPanel(type = "pills", id = "tabs",
                           ## add file preview tab
                           ##------QC plots---------
                           tabPanel("QC plots", fluidPage(
                               hr(),
                               tabsetPanel(id="QC",
                                           tabPanel(title="Violin Plots", value = "QC_panel1",
                                                    br(),
                                                    fluidRow(
                                                        column(6,
                                                               plotlyOutput("nGenePlot", width = "100%"),
                                                               plotlyOutput("mitoPlot", width = "100%"),
                                                               plotlyOutput("nUMIPlot", width = "100%")
                                                               ),
                                                        column(6,
                                                               verbatimTextOutput("name")
                                                        )
                                                    )

                                                    ),
                                           tabPanel(title="Cell and Gene Plots", value="QC_panel2",
                                                    br(),
                                                    fluidRow(
                                                      column(6,
                                                             plotlyOutput("CellPlot1", width = "50%")
                                                             ),
                                                      column(6,
                                                             plotlyOutput("CellPlot2", width = "50%")
                                                             )
                                                    )
                                           )
                               )
                           )),
                           ##------Variable Genes---------
                           tabPanel("Variable Gene Plot", fluidPage(
                               hr(),
                               textOutput("nVarGenes"),
                               fluidRow(
                                 column(4,
                                        numericInput("y.cutoff",
                                                     label = "Y Cut-off",
                                                     value = 0.5)
                                 ),
                                 column(4,
                                        numericInput("x.cutoff",
                                                     label = "X Cut-off",
                                                     value = 0.1,
                                                     min = 0)
                                 ),
                                 column(4,
                                        actionButton("findVarGenes", "Find variable genes", icon = icon("hand-pointer-o")),
                                        actionButton("doVarplot", "Plot variable genes", icon = icon("hand-pointer-o"))
                                )),
                               plotOutput("VarGenes", width = "100%")
                           )),
                           ##------PCA---------
                           tabPanel("PCA", fluidPage(
                             hr(),
                             tabsetPanel(id="Pca",
                                         tabPanel(title="PCA Plot", value="P_panel1",
                                                  br(),
                                                  fluidRow(
                                                    column(3,
                                                           actionButton("doPCA", "Run PCA", icon = icon("hand-pointer-o"))
                                                    ),
                                                    column(3,
                                                           selectInput("x.pc",
                                                                       label = "X-axis PC to use",
                                                                       choices = 1:20,
                                                                       selected = 1)
                                                           #placeholder for z-axis input
                                                    ),
                                                    column(3,
                                                           selectInput("y.pc",
                                                                       label = "Y-axis PC to use",
                                                                       choices = 1:20,
                                                                       selected = 2)
                                                    ),
                                                    column(3,
                                                           selectInput("z.pc",
                                                                       label = "Z-axis PC to use",
                                                                       choices = 1:20,
                                                                       selected = 3)
                                                    )
                                                  ),
                                                  fluidRow(
                                                      column(3 #,
                                                             ## placeholder for point size argument
                                                             #numericInput("pc.plot.size",
                                                             #             label = "Point Size:",
                                                             #             value = 1,
                                                             #             min = 0.1,
                                                             #             step = 0.5)
                                                      ),
                                                      column(3,
                                                             sliderInput("pca.plot.alpha",
                                                                         label = "Point Transparency",
                                                                         min = 0,
                                                                         max = 1,
                                                                         step = 0.1,
                                                                         value = 0.8)
                                                      ),
                                                      column(6,
                                                             uiOutput("clustUI")
                                                      )
                                                  ),
                                                  plotlyOutput("PCA2DPlot", width = "100%"),
                                                  plotlyOutput("PCA3DPlot", width = "100%")
                                         ),
                                         tabPanel(title="PC Gene Visualisation", value="P_panel2",
                                                  br(),
                                                  selectInput("select.pc",
                                                              label = "PC to plot",
                                                              choices = c(1:20)
                                                              ),
                                                  fluidRow(
                                                      column(4,
                                                             plotOutput("vizPlot", width = "100%", height = "600px")
                                                             ),
                                                      column(8,
                                                             plotOutput("PCHeatmap", width = "100%", height = "600px")
                                                      )
                                                  ),
                                                  DT::dataTableOutput("PCtable")
                                         )
                             )
                           )),
                           ##------Diagnostic---------
                           tabPanel("Diagnostic", fluidPage(
                             hr(),
                             tabsetPanel(id="diag",
                                         tabPanel(title="JackStraw", value="D_panel1",
                                                  br(),
                                                  actionButton("doJack", label = "Plot Jackstraw"),
                                                  br(),
                                                  plotOutput("Jackstraw", width = "100%")
                                         ),
                                         tabPanel(title="Elbow", value="D_panel2",
                                                  br(),
                                                  actionButton("doElbow", label = "Get Elbow Plot"),
                                                  br(),
                                                  plotOutput("Elbow", width = "100%")

                                         )
                             )
                           )),
                           ##------TSNE---------
                           tabPanel("TSNE", fluidPage(
                             hr(),
                             fluidRow(
                               column(4,
                                      numericInput("dim.used",
                                                   label = "Dimensions used",
                                                   value = 10)
                               ),
                               column(4,
                                      numericInput("max.iter",
                                                   label = "Max Iterations",
                                                   value = 2000,
                                                   min = 100)
                               ),
                               column(4,
                                      br(),
                                      actionButton("doTsne", "Run TSNE", icon = icon("hand-pointer-o")),
                                      textOutput("Tsne.done"),
                                      br()
                               )),
                             fluidRow(
                                 column(3 #,
                                        ## placeholder for point size argument
                                        #numericInput("tsne.plot.size",
                                        #             label = "Point Size:",
                                        #             value = 1,
                                        #             min = 0.1,
                                        #             step = 0.5)
                                 ),
                                 column(3,
                                        sliderInput("tsne.plot.alpha",
                                                    label = "Point Transparency",
                                                    min = 0,
                                                    max = 1,
                                                    step = 0.1,
                                                    value = 0.8)
                                 )
                             ),
                             br(),
                             fluidRow(
                                 column(10,
                                        plotlyOutput("Tsne_2d_plot", width = "100%"),
                                        plotlyOutput("Tsne_3d_plot", width = "100%")
                                 ),
                                 column(2,
                                        textOutput("selection.summary"),
                                        actionButton("create.selection", label = "Create cluster from selection"),
                                        actionButton("reset.selection", label = "Reset identities"))
                             )
                           )),
                           ##------DEGs---------
                           tabPanel("DEGs", fluidPage(
                             hr(),
                             fluidRow(
                               column(4,
                                      uiOutput("clust1")
                               ),
                               column(4,
                                      uiOutput("clust2")
                               ),
                               column(4,
                                      actionButton("doDeg", "Generate DEG plots and tables", icon = icon("hand-pointer-o"))
                               )),
                             fluidRow(
                               column(6,
                                      uiOutput("deg.gene.select"),
                                      plotlyOutput("Deg.plot", width = "100%")
                               ),
                               column(6,
                                      DT::dataTableOutput("Deg.table")
                               ))
                           ))
                           ##------END---------
               )
        )
    )
))
