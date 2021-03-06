% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{dr_scatterPlot}
\alias{dr_scatterPlot}
\title{Create interactive scatter plots of dimension reduced embeddings from Seurat object}
\usage{
dr_scatterPlot(sObj, x.axis = 1, y.axis = 2, z.axis = 3, dim = "2D",
  datatype = "pca", group.by = NULL, alpha = 0.8, point.size = 1,
  interactive = TRUE)
}
\arguments{
\item{sObj}{A Seurat object}

\item{x.axis}{An integer specifying the column number to use on the x-axis}

\item{y.axis}{An integer specifying the column number to use on the y-axis}

\item{z.axis}{An integer specifying the column number to use on the x-axis}

\item{dim}{Dimensions to plot, either "2D" or "3D"}

\item{datatype}{Which dimension reduction embeddings to use, "pca" (Default), or "tsne"}

\item{group.by}{Specifies how points should be colored}

\item{alpha}{Value between 0 and 1 specifying point opacity. Useful to set lower for large datasets.}

\item{point.size}{Numeric value to set point size. Lower settings can help prevent points being crowded out.}

\item{interactive}{If \code{TRUE}, generates plotly object, otherwise returns ggplot object.}
}
\value{
A plotly object, or ggplot object if interactive is \code{FALSE}
}
\description{
Can hover over points to see the cell ID, its current identity, and original identity.
Plots can be interactive or 3D
}
\note{
3D plots cannot be generated with interactive set to \code{FALSE}
}
\examples{

}
\author{
Matthew Myint
}
