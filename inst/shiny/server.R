## max data size
options(shiny.maxRequestSize = 1024^10)
options(shiny.launch.browser = T)

shinyServer(function(input, output, session) {
    v <- reactiveValues(scData = NULL,
                        isPCAdone = NULL,
                        isTSNEdone = NULL,
                        isClusterdone = NULL)
    celltypes <- NULL
    prePlot <- function(){
      while(names(dev.cur()) != "null device"){
        dev.off()
      }
    }
    ##-------------------Side Panel-------------------

    normMethod <- NULL

    observeEvent(input$loadButton, {
        tpmFiles <- input$tpmFiles
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
                sObj <- CreateSeuratObject(exp.data,
                              project = input$projName,
                              names.field = input$field,
                              names.delim = input$delim,
                              is.expr = input$expThres,
                              normalization.method = normMethod,
                              min.genes = input$min.genes)
                mito.genes <- grep("^MT-", rownames(sObj@data), ignore.case = TRUE, value = TRUE)
                percent.mito <- colSums(sObj@raw.data[mito.genes, ])/colSums(sObj@raw.data)
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

    ## Add option if no mito genes
    observeEvent(input$QC == "QC_panel1", {
        QC_ViolinInput <- function(){
          if(is.null(v$scData)){
            return(NULL)
          }else{
            withProgress(message="Generating Violin Plot...", value=0, {
              VlnPlot(v$scData, c("nGene", "percent.mito", "nUMI"), nCol = 1)
            })
          }
        }
        output$ViolinPlot <- renderPlot({
          QC_ViolinInput()
          }, height = 800, width = 850)
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
              pdf(filename2,
                  width=as.numeric(input$pdf_w),
                  height=as.numeric(input$pdf_h))
              print(QC_ViolinInput())
              dev.off()
            })
          }
        })
    })

    ## Cell plot
    observeEvent(input$QC == "QC_panel2", {
        if(!is.null(v$scData)){
          withProgress(message="Generating Cell/Gene Plots...", value=0, {
            gPlot1 <- GenePlot(v$scData, "nUMI", "nGene", cex.use = 1, do.hover = TRUE)
            gPlot2 <- GenePlot(v$scData, "nUMI", "percent.mito", cell.ids = WhichCells(v$scData), cex.use = 1, do.hover = TRUE)
          })
          output$CellPlot1 <- renderPlotly({
            print(gPlot1)
          })
          output$CellPlot2 <- renderPlotly({
            print(gPlot2)
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
                GenePlot(v$scData, "nUMI", "percent.mito", cell.ids = WhichCells(v$scData), cex.use = 1, do.hover = TRUE)
                dev.off()
              })
            }
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
      })
    })

    output$clustUI <- renderUI({
      if(is.null(v$isPCAdone)){
        return(NULL)
      }else{
        tagList(
          numericInput("clus.res",
                       label = "Cluster Resolution",
                       value = 0.6,
                       min = 0.1,
                       step = 0.1),
          br(),
          actionButton("findCluster", "Find Clusters", icon = icon("hand-pointer-o")),
          textOutput("cluster.done")
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

    output$pca_plotspace <- renderUI({
      if(is.null(v$isPCAdone)){
        return(NULL)
      }else{
        plotlyOutput("PCAPlot", width = "100%")
      }
    })

    observeEvent(input$Pca == "P_panel1", {
      pcaPlotInput <- function(){
        if(is.null(v$scData) || !v$isPCAdone){
          return(NULL)
        }else{
          withProgress(message="Generating PCA Plot...", value=0, {
            PCAPlot(v$scData, dim.1 = input$x.pc, dim.2 = input$y.pc, pt.size = 2, do.hover = TRUE)
          })
        }
      }
      output$PCAPlot <- renderPlotly({
        pcaPlotInput()
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
            prePlot()
            pdf(filename2,
                width=as.numeric(input$pdf_w),
                height=as.numeric(input$pdf_h))
            PCAPlot(v$scData, dim.1 = input$x.pc, dim.2 = input$y.pc, pt.size = 2)
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
    })

    # Viz plot
    observeEvent(input$Pca == "P_panel2", {
      if(!is.null(v$scData)){
        v$scData <- ProjectPCA(v$scData)
      }
      VizInput <- function(){
        if(is.null(v$scData)){
          return(NULL)
        }else{
          withProgress(message="Generating Viz Plot...", value=0, {
            VizPCA(v$scData)
          })
        }
      }
      output$vizPlot <- renderPlot({
        VizInput()
      }, height = 800, width = 850)
      observeEvent(input$PDFe, {
        if(!is.null(v$scData)){
          withProgress(message="Downloading plot PDF files...", value=0, {
            print(getwd())
            pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
            if(!dir.exists(pdfDir)){
              dir.create(pdfDir)
            }
            filename2 <- paste0(pdfDir, .Platform$file.sep,"Viz_plot_", Sys.Date(), ".pdf")
            i = 0
            while(file.exists(filename2)){
              filename2 <- paste0(pdfDir, .Platform$file.sep,
                                  "Viz_plot_",
                                  Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
              i = i + 1;
            }
            prePlot()
            pdf(filename2,
                width=as.numeric(input$pdf_w),
                height=as.numeric(input$pdf_h))
            VizInput()
            dev.off()
          })
        }
      })
    })

    # PC heatmap
    observeEvent(input$Pca == "P_panel3", {
      PCHmInput <- function(){
        if(is.null(v$scData)){
          return(NULL)
        }else{
          withProgress(message="Generating PC Heatmap...", value=0, {
            PCHeatmap(v$scData, pc.use = input$PCused, do.balanced = TRUE, label.columns = TRUE, use.full = FALSE, srtCol = 45, offsetRow = -0.5, cexCol = 0.5, offsetCol = -0.5, key = FALSE)
          })
        }
      }
      output$PCHeatmap <- renderPlot({
        PCHmInput()
      }, height = 800, width = 850)
      observeEvent(input$PDFf, {
        if(!is.null(v$scData)){
          withProgress(message="Downloading plot PDF files...", value=0, {
            print(getwd())
            pdfDir <- paste0(getwd(), .Platform$file.sep, "Seurat_results/Generated_reports_", Sys.Date())
            if(!dir.exists(pdfDir)){
              dir.create(pdfDir)
            }
            filename2 <- paste0(pdfDir, .Platform$file.sep,"PC_Heatmap_PC", input$PCused, "_", Sys.Date(), ".pdf")
            i = 0
            while(file.exists(filename2)){
              filename2 <- paste0(pdfDir, .Platform$file.sep,"PC_Heatmap_PC", input$PCused, "_", Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
              i = i + 1;
            }
            prePlot()
            pdf(filename2,
                width=as.numeric(input$pdf_w),
                height=as.numeric(input$pdf_h))
            PCHmInput()
            dev.off()
          })
        }
      })
    })

    ##---------------Significant PCs tabset-------------------

    # Jackstraw
    observeEvent(input$doJack, { #put under control of run button again
      JackInput <- function(){
        if(is.null(v$scData)){
          return(NULL)
        }else{
          withProgress(message="Generating Jackstraw Plot...", value=0, {
            v$scData <- JackStraw(v$scData, num.replicate = 100, do.print = TRUE)
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
    output$perplex.option <- renderUI({
        if(is.null(v$isPCAdone)){
            return(NULL)
        }else{
            ##perplexity test
            n.cells <- isolate(nrow(v$scData@dr$pca@cell.embeddings))
            max.perplex <- as.integer((n.cells - 1)/3)
            numericInput("perplexity",
                         label = "Perplexity",
                         value = if(max.perplex <30) max.perplex else 30,
                         min = 0,
                         max = max.perplex)
        }
    })

    observeEvent(input$doTsne, {
      withProgress(message = "Running tSNE...", value = 0.3, {
        v$scData <- RunTSNE(v$scData, dims.use = 1:input$dim.used, max_iter = input$max.iter,
                            do.fast = TRUE, perplexity = input$perplexity)
        output$Tsne.done <- renderText(paste0("TSNE done!"))
        v$isTSNEdone <- TRUE
      })
    })

    observeEvent(input$doTsnePlot, {
      TsneInput <- function(){
        if(is.null(v$scData)){
          return(NULL)
        }else{
          withProgress(message="Generating TSNE Plot...", value=0, {
            TSNEPlot(v$scData, pt.size = 2, do.hover = TRUE)
          })
        }
      }
      output$Tsne.plot <- renderPlotly({
        TsneInput()
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
            prePlot()
            pdf(filename2,
                width=as.numeric(input$pdf_w),
                height=as.numeric(input$pdf_h))
            TSNEPlot(v$scData, pt.size = 2)
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
    })

    ##---------------DEGs tabset-------------------
    output$clust1 <- renderUI({
      if(is.null(v$scData)){
        return(NULL)
      }else{
        celltypes <- levels(v$scData@ident)
        selectInput('c1', 'Choose cluster of interest:',
                    choices = celltypes,
                    selected = celltypes[1],
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
          vlnMarkers <- rownames(ips.markers)
          if(length(vlnMarkers) < 5 ){
            vlnMarkers <- vlnMarkers[!is.na(vlnMarkers)]
          }else{
            vlnMarkers <- vlnMarkers[1:5]
          }
        })
        DegPlotInput <- function(){
          if(is.null(v$scData)){
            return(NULL)
          }else{
            withProgress(message="Generating DEG Plots...", value=0, {
              VlnPlot(v$scData, vlnMarkers, nCol = 1)
            })
          }
        }
        output$Deg.plot <- renderPlot({
          DegPlotInput()
        }, height = 800, width = 400)
        output$Deg.table <- renderTable(print(head(ips.markers, 15)), rownames = TRUE, digits = -1)

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
              prePlot()
              pdf(filename2,
                  width=as.numeric(input$pdf_w),
                  height=as.numeric(input$pdf_h))
              print(DegPlotInput())
              dev.off()
              write.csv(ips.markers, file = paste0(pdfDir, .Platform$file.sep,"DEG_plot_", input$c1, "vs", input$c2, "_", Sys.Date(), ".csv"))
            })
          }
        })
      }
    })
    ##---------------Summary tab

    ##------Clean up when ending session----
    session$onSessionEnded(function(){
      prePlot()
    })
})


