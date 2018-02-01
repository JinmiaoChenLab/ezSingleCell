require(shiny)
require(Seurat)
require(dplyr)
require(Matrix)
require(plotly)

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
                                        ".fcs")),
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
                                         value = 2,
                                         min = 1)
                     )
                   ),
                   numericInput(inputId = "expThres",
                                label = "Expression cut off",
                                value = 1,
                                min = 0.1),
                   numericInput(inputId = "min.genes",
                                label = "Minimum genes",
                                value = 200,
                                min = 1),
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
                                  actionButton("PDFe", "Download Viz Plot in PDF", icon = icon("download"))
                 ),
                 conditionalPanel(" input.Pca == 'P_panel3' && input.tabs == 'PCA' ",
                                  actionButton("PDFf", "Download PCA Heatmap in PDF", icon = icon("download"))
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
                                                    plotOutput("ViolinPlot", width = "100%")
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
                                         tabPanel(title="Plot", value="P_panel1",
                                                  br(),
                                                  fluidRow(
                                                    column(2,
                                                           actionButton("doPCA", "Run PCA", icon = icon("hand-pointer-o"))
                                                    ),
                                                    column(5,
                                                           numericInput("x.pc",
                                                                        label = "X-axis PC to use",
                                                                        value = 1),
                                                           numericInput("y.pc",
                                                                        label = "Y-axis PC to use",
                                                                        value = 2)
                                                    ),
                                                    column(5,
                                                           uiOutput("clustUI")

                                                    )),
                                                  uiOutput("pca_plotspace")
                                         ),
                                         tabPanel(title="Viz", value="P_panel2",
                                                  br(),
                                                  plotOutput("vizPlot", width = "100%")
                                         ),
                                         tabPanel(title="PC Heatmap", value="P_panel3",
                                                  br(),
                                                  fluidRow(
                                                    column(2,
                                                           numericInput("PCused",
                                                                        label = "PC to select for heatmap",
                                                                        value = 10,
                                                                        min = 1)
                                                    ),
                                                    column(10,
                                                           plotOutput("PCHeatmap", width = "100%")
                                                    ))
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
                               column(3,
                                      numericInput("dim.used",
                                                   label = "Dimensions used",
                                                   value = 10)
                               ),
                               column(3,
                                      numericInput("max.iter",
                                                   label = "Max Iterations",
                                                   value = 2000,
                                                   min = 100)
                               ),
                               column(3,
                                      uiOutput("perplex.option")
                               ),
                               column(3,
                                      br(),
                                      actionButton("doTsne", "Run TSNE", icon = icon("hand-pointer-o")),
                                      textOutput("Tsne.done"),
                                      br(),
                                      actionButton("doTsnePlot", "Update TSNE plot", icon = icon("hand-pointer-o"))

                               )),
                             plotlyOutput("Tsne.plot", width = "100%")
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
                                      plotOutput("Deg.plot", width = "50%")
                               ),
                               column(6,
                                      tableOutput("Deg.table")
                               ))
                           ))
                           ##------END---------
               )
        )
    )
))
