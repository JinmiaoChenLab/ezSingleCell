#' Create interactive violin plot from Seurat object
#'
#' Can hover over points to see the cell ID, its current identity, original identity, and the value of the parameter being plotted.
#' @return A plotly object
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @param sObj The seurat object
#' @param type The parameter to be visualised. Currently applies to metadata (nGene, nUmi, percent.mito)
#' @param interactive If \code{TRUE}, output will be plotly object, else a ggplot object.
#' @author Matthew Myint
#' @export
#' @examples
#'
violin_plot <- function(sObj, type, interactive = TRUE, event.data = NULL) {
    metadata <- sObj@meta.data
    y.val <- switch(type,
                    "nGene" = metadata$nGene,
                    "nUMI" = metadata$nUMI,
                    #type for gene values,
                    "percent.mito" = metadata$percent.mito)
    if(type %in% rownames(sObj@raw.data)){
        y.val <- as.vector(FetchData(sObj, type))
        metadata <- cbind(metadata, y.val)
    }
    orig.ident <- factor(metadata$orig.ident)
    current.ident <- factor(sObj@ident)
    keys <- row.names(metadata)
    set.seed(2017)
    p <-
        ggplot(
            metadata,
            aes(
                x = current.ident,
                y = y.val
            )
        ) +
        geom_violin(trim = T,
                    scale = "count",
                    width = 0.7,
                    aes(fill = current.ident)) +
        ggtitle(type) +
        ## using log scale for y axis
        # scale_y_log10() +
        geom_jitter(shape = 21,
                    size = 1,
                    cex = 0.5,
                    colour = "black",
                    aes(text = keys,
                        key = keys,
                        fill = orig.ident)
                    ,
                    width = 0.3) +
        theme(legend.position = "none",
              axis.title = element_blank())
    pp <- ggplotly(p, tooltip = c("text", "x", "fill", "y")) %>%
        layout(xaxis = list(title = NULL),
               yaxis = list(title = NULL),
               dragmode = "select")
    if(!is.null(event.data)){
        m <- metadata[keys %in% event.data,]
        pp <- add_markers(pp, data = m, alpha = 0.3, color = I("yellow"), size = 1.5)
    }
    if(interactive){
        return(pp)
    }else{
        return(p)
    }
}

#' Create interactive scatter plots of dimension reduced embeddings from Seurat object
#'
#' Can hover over points to see the cell ID, its current identity, and original identity.
#' Plots can be interactive or 3D
#' @param sObj A Seurat object
#' @param x.axis An integer specifying the column number to use on the x-axis
#' @param y.axis An integer specifying the column number to use on the y-axis
#' @param z.axis An integer specifying the column number to use on the x-axis
#' @param dim Dimensions to plot, either "2D" or "3D"
#' @param datatype Which dimension reduction embeddings to use, "pca" (Default), or "tsne"
#' @param group.by Specifies how points should be colored
#' @param alpha Value between 0 and 1 specifying point opacity. Useful to set lower for large datasets.
#' @param interactive If \code{TRUE}, generates plotly object, otherwise returns ggplot object
#' @return A plotly object, or ggplot object if interactive is \code{FALSE}
#' @note 3D plots cannot be generated with interactive set to \code{FALSE}
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @importFrom scales hue_pal
#' @author Matthew Myint
#' @export
#' @examples
#'
dr_scatterPlot <- function(sObj, x.axis = 1, y.axis = 2, z.axis = 3,
                           dim = "2D", datatype = "pca", group.by = NULL,
                           alpha = 0.8, interactive = TRUE){
    if(dim == "3D" && !interactive){
        warning("3D plot can only be interactive, setting interactive to TRUE!")
        interactive <- TRUE
    }
    data <- switch(datatype,
                   "pca" = sObj@dr$pca@cell.embeddings,
                   "tsne" = sObj@dr$tsne@cell.embeddings)
    if(is.null(group.by)){
        group.by <- sObj@ident
    }
    Cluster <- factor(group.by)
    data <- cbind(as.data.frame(data), group.by)
    key <- rownames(data)
    if(dim == "2D"){
        Coordinates <- paste0("(", round(data[, x.axis], 2), ",", round(data[, y.axis], 2), ")")
        p <- ggplot(data, aes(x = data[, x.axis],
                           y = data[, y.axis],
                           label = Coordinates)) +
            geom_point(aes(colour = Cluster,
                           text = key,
                           key = key),
                       alpha = alpha) +
            theme(legend.position = "none",
                  axis.title.y = element_text(angle = 90))+
            labs(x = colnames(data)[x.axis],
                 y = colnames(data)[y.axis])

        pp <- ggplotly(p, tooltip = c("text", "colour", "label")) %>%
            layout(scene = list(
                dragmode = "select"
            ))
        if(interactive){
            return(pp)
        }else{
            return(p)
        }
    }else if(dim == "3D"){
      Coordinates <- paste0("(", round(data[, x.axis], 2), ",",
                            round(data[, y.axis], 2), ",",
                            round(data[, z.axis], 2), ")")
      pp <- plot_ly(
            data,
            x = data[, x.axis],
            y = data[, y.axis],
            z = data[, z.axis],
            color = Cluster,
            colors = hue_pal()(length(levels(group.by))),
            key = key,
            alpha = alpha,
            mode = "markers",
            type = "scatter3d",
            hoverinfo = "text",
            text = ~paste(key,
                          "<br> Cluster: ", Cluster,
                          "<br>", Coordinates)
        )  %>%
            layout(scene = list(
                xaxis = list(title = colnames(data)[x.axis]),
                yaxis = list(title = colnames(data)[y.axis]),
                zaxis = list(title = colnames(data)[z.axis]),
                dragmode = "select"
            ))
    }
    return(pp)
}
