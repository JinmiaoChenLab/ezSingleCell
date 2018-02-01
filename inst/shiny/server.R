## max data size
options(shiny.maxRequestSize = 1024^10)
options(shiny.launch.browser = T)

shinyServer(function(input, output, session) {
    v <- reactiveValues(scData = NULL,
                        isPCAdone = NULL,
                        isTSNEdone = NULL,
                        isClusterdone = NULL,
                        pcGenes = NULL,
                        plotlySelection = NULL,
                        ips.markers = NULL)
    celltypes <- NULL
    prePlot <- function(){
      while(names(dev.cur()) != "null device"){
        dev.off()
      }
    }
    observe({
        #s <- event_data("plotly_selected")
        #cells <- s[["key"]]
        v$plotlySelection <- event_data("plotly_selected")[["key"]]
    })
    ##-------------------Side Panel-------------------

    normMethod <- NULL

    observeEvent(input$loadButton, {
        tpmFiles <- input$tpmFiles
        annoFile <- input$cellAnnoFiles
        if(input$norm){
          normMethod <- "LogNormalize"
        }
        if (is.null(tpmFiles)){
            v$scData <- NULL
        }else{
            withProgress(message="Loading and Processing Data...", value=0, {
                print(tpmFiles$datapath)
                print(tpmFiles$name)
                print(file.exists(paste(tpmFiles$datapath[1], "/", tpmFiles$name[1], sep="")))
                exp.data <- read.table(tpmFiles$datapath,
                                       sep="\t", header=TRUE, row.names=1, stringsAsFactors = FALSE)
                anno.data <- read.table(annoFile$datapath[1], header = T, sep = "\t",
                                        row.names = 1, stringsAsFactors = FALSE)
                anno.data$combined <- paste(anno.data, anno.data, anno.data, sep = "_")
                incProgress(0.5, "Creating Seurat Object")
                sObj <- CreateSeuratObject(exp.data,
                              project = input$projName,
                              names.field = input$field,
                              names.delim = input$delim,
                              is.expr = input$expThres,
                              normalization.method = normMethod,
                              min.genes = input$min.genes,
                              min.cells = input$min.cells)
                mito.genes <- grep("^MT-", rownames(sObj@data), ignore.case = TRUE, value = TRUE)
                percent.mito <- colSums(sObj@raw.data[mito.genes, ])/colSums(sObj@raw.data)
                incProgress(0.5, "Adding metadata")
                sObj <- AddMetaData(sObj, percent.mito, "percent.mito")
                v$scData <- sObj
            })
        }
        dir.create("Seurat_results")
    })

    observeEvent(input$reset, {
      session$reload()
      print("Reset done")
    })

    observeEvent(input$saveButton, {
        if(!is.null(input$tpmFiles)){
            withProgress(message="Saving Results...", value=0, {
                print(getwd())
                dir.create("Seurat_results")
                resultDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results")
                filename <- paste0(resultDir, .Platform$file.sep, v$scData@project.name, "_", Sys.Date())
                sObj <- v$scData
                save(sObj, file= paste0(resultDir, .Platform$file.sep, sObj@project.name, "_", Sys.Date(), ".Robj"))
            })
          ## open the results directory
          opendir(resultDir)
        }
    })

    output$logo <- renderImage({
      return(list(
        src = "inst/extdata/logo.png",
        contentType = "image/png",
        alt = "Singapore Immunology Network"
      ))
    }, deleteFile = FALSE)

    opendir <- function(dir = getwd()){
      if (.Platform['OS.type'] == "windows"){
        shell.exec(dir)
      } else {
        system(paste(Sys.getenv("R_BROWSER"), dir))
      }
    }

    observeEvent(input$OpenDir, {
      resultDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results")
      if(!dir.exists(resultDir)){
        dir.create("Seurat_results")
      }
      pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
      if(dir.exists(pdfDir)){
        opendir(pdfDir)
      }else{
        warning("No reports created yet!")
        dir.create(pdfDir)
      }
    })

    ##---------------QC tabset-------------------

    output$nGenePlot <- renderPlotly({
        if(is.null(v$scData)){
            plotly_empty()
        }else{
            violin_plot(v$scData, "nGene")
        }
    })

    output$mitoPlot <- renderPlotly({
        if(is.null(v$scData)){
            plotly_empty()
        }else{
            violin_plot(v$scData, "percent.mito")
        }
    })

    output$nUMIPlot <- renderPlotly({
        if(is.null(v$scData)){
            plotly_empty()
        }else{
            violin_plot(v$scData, "nUMI")
        }
    })

    output$name <- renderPrint({
        s <- event_data("plotly_selected")
        c(s[["key"]], class(s[["key"]]))
    })

    observeEvent(input$PDFa, {
        if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"QC_violin_plot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,
                                        "QC_violin_plot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                prePlot()
                nG <- violin_plot(v$scData, "nGene", interactive = FALSE)
                pM <- violin_plot(v$scData, "percent.mito", interactive = FALSE)
                nU <- violin_plot(v$scData, "nUMI", interactive = FALSE)
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                print(nG)
                print(pM)
                print(nU)
                dev.off()
            })
        }
    })

    ## Cell plot

    output$CellPlot1 <- renderPlotly({
        if(is.null(v$scData)){
            plotly_empty()
        }else{
            print(GenePlot(v$scData, "nUMI", "nGene", cex.use = 1, do.hover = TRUE))
        }
    })

    output$CellPlot2 <- renderPlotly({
        if(is.null(v$scData)){
            plotly_empty()
        }else{
            print(GenePlot(v$scData, "nUMI", "percent.mito", cell.ids = WhichCells(v$scData), cex.use = 1, do.hover = TRUE))
        }
    })

    observeEvent(input$PDFb, {
        if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"QC_cell_plot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,
                                        "QC_cell_plot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                prePlot()
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                GenePlot(v$scData, "nUMI", "nGene", cex.use = 1)
                GenePlot(v$scData, "nUMI", "percent.mito", cell.ids = WhichCells(v$scData), cex.use = 1)
                dev.off()
            })
        }
    })

    ##---------------Variable Genes tabset-------------------
    observeEvent(input$findVarGenes, {
      withProgress(message = "Finding variable genes...", value = 0, {
        v$scData <- FindVariableGenes(v$scData,
                                      mean.function = ExpMean,
                                      dispersion.function = LogVMR,
                                      x.low.cutoff = input$x.cutoff,
                                      y.cutoff = input$y.cutoff,
                                      do.plot = FALSE)
        incProgress(0.5)
        VarGeneText <- paste0("Number of variable genes: ", length(v$scData@var.genes))
        output$nVarGenes <- renderText(VarGeneText)
      })
    observeEvent(input$doVarplot,{
      varGenePlotInput <- function(){
        if(is.null(v$scData)){
            return(NULL)
        }else{
            withProgress(message="Plotting variable genes...", value=0, {
              VariableGenePlot(v$scData,
                               x.low.cutoff = input$x.cutoff,
                               y.cutoff = input$y.cutoff,
                               do.contour = FALSE)
            })
        }
      }

    output$VarGenes <- renderPlot({
          varGenePlotInput()
        }, height = 800, width = 850)
    observeEvent(input$PDFc, {
          if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
              print(getwd())
              pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
              if(!dir.exists(pdfDir)){
                dir.create(pdfDir)
              }
              filename2 <- paste0(pdfDir, .Platform$file.sep,"Var_genes_plot_", Sys.Date(), ".pdf")
              i = 0
              while(file.exists(filename2)){
                filename2 <- paste0(pdfDir, .Platform$file.sep,
                                    "Var_genes_plot_",
                                    Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                i = i + 1;
              }
              prePlot()
              pdf(filename2,
                  width=as.numeric(input$pdf_w),
                  height=as.numeric(input$pdf_h))
              varGenePlotInput()
              mtext(VarGeneText)
              dev.off()
              txtfile <- sub("Var_genes_plot_", "Var_gene_list_", filename2)
              txtfile <- sub(".pdf", ".txt", txtfile)
              write(v$scData@var.genes, file = txtfile)
            })
          }
        })
      })
    })

    ##---------------PCA tabset-------------------
    # PCA plot
    observeEvent(input$doPCA, {
      withProgress(message = "Scaling Data...", value = 0,{
        v$scData <- ScaleData(v$scData)
        incProgress(0.5, message = "Running PCA...")
        v$scData <- RunPCA(v$scData, pc.genes = v$scData@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
        v$isPCAdone <- TRUE
        v$scData <- ProjectPCA(v$scData)
        incProgress(0.4, message = "Getting list of PC genes...")
        pc.table <- list()
        for(i in 1:20){
            pcg <- DimTopGenes(v$scData, dim.use = i, do.balanced = TRUE, num.genes = 30)
            pc.table[[i]] <- pcg
        }
        pc.table <- as.data.frame(pc.table, col.names = paste0("PC", 1:20))
        v$pcGenes <- pc.table
      })
    })

    output$clustUI <- renderUI({
      if(is.null(v$isPCAdone)){
        return(NULL)
      }else{
        tagList(
          fluidRow(
              column(6,
                     numericInput("clus.res",
                                  label = "Cluster Resolution",
                                  value = 0.6,
                                  min = 0.1,
                                  step = 0.1)
                     ),
              column(6,
                     actionButton("findCluster", "Find Clusters", icon = icon("hand-pointer-o")),
                     textOutput("cluster.done")
                     )
          )
        )
      }
    })

    observeEvent(input$findCluster, {
      withProgress(message = "Finding clusters...", value = 0.3, {
        v$scData <- FindClusters(v$scData, reduction.type = "pca", dims.use = 1:input$dim.used,
                                 resolution = input$clus.res, print.output = 0, save.SNN = TRUE)
        output$cluster.done <- renderText(paste0("Clustering done!"))
        v$isClusterdone <- TRUE
      })
    })

    output$PCA2DPlot <- renderPlotly({
        if(is.null(v$isPCAdone)){
            plotly_empty()
        }else{
            withProgress(message="Generating 2D PCA Plot...", value=0, {
                dr_scatterPlot(v$scData, x.axis = as.numeric(input$x.pc),
                               y.axis = as.numeric(input$y.pc),
                               alpha = input$pca.plot.alpha,
                               dim = "2D", datatype = "pca")
            })
        }
    })

    output$PCA3DPlot <- renderPlotly({
        if(is.null(v$isPCAdone)){
            plotly_empty()
        }else{
            withProgress(message="Generating 3D PCA Plot...", value=0, {
                dr_scatterPlot(v$scData, x.axis = as.numeric(input$x.pc),
                               y.axis = as.numeric(input$y.pc),
                               z.axis = as.numeric(input$z.pc),
                               alpha = input$pca.plot.alpha,
                               dim = "3D", datatype = "pca")
            })
        }
    })

    observeEvent(input$PDFd, {
        if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"PCA_plot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,
                                        "PCA_plot_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                pcaplot <- dr_scatterPlot(v$scData, x.axis = as.numeric(input$x.pc),
                                          y.axis = as.numeric(input$y.pc),
                                          alpha = input$pca.plot.alpha,
                                          dim = "2D", datatype = "pca", interactive = FALSE)
                prePlot()
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                print(pcaplot)
                dev.off()
            })
            withProgress(message="Downloading PCA coordinates...", value=0.5, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"pca_", Sys.Date(), ".txt")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,
                                        "pca_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".txt");
                    i = i + 1;
                }
                write.table(v$scData@dr$pca@cell.embeddings, file = filename2)
            })
            withProgress(message="Downloading cluster IDs...", value=0.9, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"cluster_", Sys.Date(), ".txt")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"cluster_", Sys.Date(), "_", sprintf("%03d", i + 1), ".txt");
                    i = i + 1;
                }
                write.table(v$scData@ident, file = filename2)
            })
        }
    })

    # Viz plot

    output$vizPlot <- renderPlot({
        if(is.null(v$scData) && v$isPCAdone){
            return(NULL)
        }else{
            VizPCA(v$scData, pcs.use = as.numeric(input$select.pc))
        }
    }, width = 600)

    output$PCHeatmap <- renderPlot({
        if(is.null(v$scData) && v$isPCAdone){
            return(NULL)
        }else{
            PCHeatmap(v$scData, pc.use = as.numeric(input$select.pc), do.balanced = TRUE, label.columns = TRUE, use.full = FALSE, srtCol = 45, offsetRow = -0.5, cexCol = 0.5, offsetCol = -0.5, key = FALSE)
        }
    })

    output$PCtable <- DT::renderDataTable({
        if(is.null(v$scData) && v$isPCAdone){
            return(NULL)
        }else{
            v$pcGenes
        }
    }, options = list(scrollX = TRUE))

    observeEvent(input$PDFe, {
        if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"Viz_Heatmap_plots_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,
                                        "Viz_Heatmap_plots_",
                                        Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                prePlot()
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                isolate({
                    VizPCA(v$scData, pcs.use = as.numeric(input$select.pc))
                    PCHeatmap(v$scData, pc.use = as.numeric(input$select.pc), do.balanced = TRUE, label.columns = TRUE, use.full = FALSE, srtCol = 45, offsetRow = -0.5, cexCol = 0.5, offsetCol = -0.5, key = FALSE)
                })
                dev.off()
                pcGenes <- v$pcGenes
                write.csv(v$pcGenes, file = paste0(pdfDir, .Platform$file.sep,"PC_genes_", Sys.Date(), ".csv"))
            })
        }
    })


    ##---------------Significant PCs tabset-------------------

    # Jackstraw
    observeEvent(input$doJack, {
      JackInput <- function(){
        if(is.null(v$scData)){
          return(NULL)
        }else{
          withProgress(message="Running Jackstraw...", value=0, {
            v$scData <- JackStraw(v$scData, num.replicate = 100, do.print = TRUE)
            incProgress(0.7, "Plotting Jackstraw...")
            JackStrawPlot(v$scData, PCs = 1:12)
          })
        }
      }
      output$Jackstraw <- renderPlot({
        JackInput()
      }, height = 800, width = 850)
      observeEvent(input$PDFg, {
        if(!is.null(v$scData)){
          withProgress(message="Downloading plot PDF files...", value=0, {
            print(getwd())
            pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
            if(!dir.exists(pdfDir)){
              dir.create(pdfDir)
            }
            filename2 <- paste0(pdfDir, .Platform$file.sep,"Jackstraw_plot_", Sys.Date(), ".pdf")
            i = 0
            while(file.exists(filename2)){
              filename2 <- paste0(pdfDir, .Platform$file.sep,"Jackstraw_plot_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
              i = i + 1;
            }
            prePlot()
            pdf(filename2,
                width=as.numeric(input$pdf_w),
                height=as.numeric(input$pdf_h))
            JackStrawPlot(v$scData, PCs = 1:12)
            dev.off()
          })
        }
      })
    })

    # Elbow
    observeEvent(input$doElbow, {
      ElbowInput <- function(){
        if(is.null(v$scData)){
          return(NULL)
        }else{
          withProgress(message="Generating Elbow Plot...", value=0, {
            PCElbowPlot(v$scData)
          })
        }
      }
      output$Elbow <- renderPlot({
        ElbowInput()
      }, height = 800, width = 850)
      observeEvent(input$PDFh, {
        if(!is.null(v$scData)){
          withProgress(message="Downloading plot PDF files...", value=0, {
            print(getwd())
            pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
            if(!dir.exists(pdfDir)){
              dir.create(pdfDir)
            }
            filename2 <- paste0(pdfDir, .Platform$file.sep,"Elbow_plot_", Sys.Date(), ".pdf")
            i = 0
            while(file.exists(filename2)){
              filename2 <- paste0(pdfDir, .Platform$file.sep,"Elbow_plot_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
              i = i + 1;
            }
            prePlot()
            pdf(filename2,
                width=as.numeric(input$pdf_w),
                height=as.numeric(input$pdf_h))
            print(ElbowInput())
            dev.off()
          })
        }
      })
    })

    ##---------------TSNE tabset-------------------
    observeEvent(input$doTsne, {
      withProgress(message = "Running tSNE...", value = 0.3, {
        v$scData <- RunTSNE(v$scData, dims.use = 1:input$dim.used, max_iter = input$max.iter, do.fast = TRUE, dim.embed = 3)
        output$Tsne.done <- renderText(paste0("TSNE done!"))
        v$isTSNEdone <- TRUE
      })
    })

    output$Tsne_2d_plot <- renderPlotly({
        if(is.null(v$scData) || is.null(v$isTSNEdone)){
            plotly_empty()
        }else{
            withProgress(message="Generating TSNE 2D Plot...", value=0, {
                dr_scatterPlot(v$scData, x.axis = 1, y.axis = 2, dim = "2D",
                               datatype = "tsne", alpha = input$tsne.plot.alpha, interactive = TRUE)
            })
        }
    })

    output$Tsne_3d_plot <- renderPlotly({
        if(is.null(v$scData) || is.null(v$isTSNEdone)){
            plotly_empty()
        }else{
            withProgress(message="Generating TSNE 3D Plot...", value=0, {
                dr_scatterPlot(v$scData, x.axis = 1, y.axis = 2, z.axis = 3, dim = "3D",
                               datatype = "tsne", alpha = input$tsne.plot.alpha, interactive = TRUE)
            })
        }
    })

    output$selection.summary <- renderText({
        if(is.null(v$plotlySelection)){
            return(NULL)
        }else{
            t <- paste0(length(v$plotlySelection), " cells selected")
            t
        }
    })

    observeEvent(input$create.selection, {
        ## stash old identity
        v$scData <- StashIdent(object = v$scData, save.name = 'cluster.ident')
        v$scData <- SetIdent(object = v$scData,
                             cells.use = v$plotlySelection,
                             ident.use = 'Selection'
        )
        updateTabsetPanel(session, "tabs", selected = "DEGs")
    })

    observeEvent(input$reset.selection, {
        v$scData <- SetAllIdent(object = v$scData,  id = 'cluster.ident')
        #event_data("plotly_select") <- NULL
        v$plotlySelection <- NULL
    })

    observeEvent(input$PDFi, {
        if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"TSNE_plot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"TSNE_plot_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                tsneplot <- dr_scatterPlot(v$scData, x.axis = 1, y.axis = 2, dim = "2D",
                                           datatype = "tsne", alpha = 0.8, interactive = TRUE)
                prePlot()
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                print(tsneplot)
                dev.off()
            })
            withProgress(message="Downloading tSNE coordinates...", value=0.6, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"tsne_", Sys.Date(), ".txt")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"tsne_", Sys.Date(), "_", sprintf("%03d", i + 1), ".txt");
                    i = i + 1;
                }
                write.table(v$scData@dr$tsne@cell.embeddings, file = filename2)
            })
        }
    })

    ##---------------DEGs tabset-------------------
    output$clust1 <- renderUI({
      if(is.null(v$scData)){
        return(NULL)
      }else{
        celltypes <- levels(v$scData@ident)
        selectInput('c1', 'Choose cluster of interest:',
                    choices = celltypes,
                    selected = if("Selection" %in% celltypes) "Selection" else celltypes[1],
                    selectize = FALSE,
                    width = "100%")
      }
    })
    output$clust2 <- renderUI({
      if(is.null(v$scData)){
        return(NULL)
      }else{
        celltypes <- levels(v$scData@ident)
        selectInput('c2', 'Choose cluster to compare to:',
                    choices = c("All", setdiff(celltypes, input$c1)),
                    selected = "All",
                    selectize = FALSE,
                    width = "100%")
      }
    })

    observeEvent(input$doDeg, {
      if(is.null(v$scData)){
        return(NULL)
      }else{
        withProgress(message="Finding DEGs...", value=0, {
          if(input$c2 == "All"){
            ips.markers <- FindMarkers(v$scData, ident.1 = input$c1, thresh.use = 2)
          }else{
            ips.markers <- FindMarkers(v$scData, ident.1 = input$c1, ident.2 = input$c2, thresh.use = 2)
          }
          ips.markers$adj_p_val <- p.adjust(ips.markers$p_val, method = "BH")
          v$ips.markers <- ips.markers
        })
      }
    })

    output$deg.gene.select <- renderUI({
        if(is.null(v$ips.markers)){
            return(NULL)
        }else{
            selectInput("deg.gene", label = "Gene to visualise",
                        choices = rownames(v$ips.markers))
        }
    })

    output$Deg.plot <- renderPlotly({
        if(is.null(v$ips.markers)){
            return(NULL)
        }else{
            withProgress(message="Generating DEG Plot...", value=0, {
                violin_plot(v$scData, input$deg.gene)
            })
        }
    })

    output$Deg.table <- DT::renderDataTable({
        if(is.null(v$scData) && v$isPCAdone){
            return(NULL)
        }else{
            signif(v$ips.markers, 4)
        }
    }, options = list(scrollX = TRUE, scrollY = "400px"))

    observeEvent(input$PDFj, {
        if(!is.null(v$scData)){
            withProgress(message="Downloading plot PDF files...", value=0, {
                print(getwd())
                pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
                if(!dir.exists(pdfDir)){
                    dir.create(pdfDir)
                }
                filename2 <- paste0(pdfDir, .Platform$file.sep,"DEG_plot_", input$c1, "vs", input$c2, "_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename2)){
                    filename2 <- paste0(pdfDir, .Platform$file.sep,"DEG_plot_", input$c1, "vs", input$c2, "_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                    i = i + 1;
                }
                degVln <- violin_plot(v$scData, input$deg.gene, interactive = FALSE)
                prePlot()
                pdf(filename2,
                    width=as.numeric(input$pdf_w),
                    height=as.numeric(input$pdf_h))
                print(degVln)
                dev.off()
                write.csv(v$ips.markers, file = paste0(pdfDir, .Platform$file.sep,"DEG_table_", input$c1, "vs", input$c2, "_", Sys.Date(), ".csv"))
            })
        }
    })
    ##---------------Summary tab

    ##------Clean up when ending session----
    session$onSessionEnded(function(){
      prePlot()
    })
})


