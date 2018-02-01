\name{NEWS}
\title{News for Package \pkg{ezSingleCell}}

\section{To do}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item More WhichCells implementation
      \item Add options to specify DEG cut off
      \item Get DEG by adj. p-val, sort by avg diff
      \item Add option for TSNE complexity
      \item Names field to show examples, not just number (render UI)
      \item Allow cluster colour to change depending on identity
      \item TPM preview/editing (preprocessing)
      \item Guide on acceptable TPM table formats
      \item UI overhaul utilising shinyDashboard
      \item Allow highlighting cells in interactive plots for DEG analysis
    }
  }
}

\section{Version 0.0.2 (2018-02-01)}{
  \subsection{Fixes}{
    \itemize{
      \item Fixed plots downloading in triplicate
      \item DESCRIPTION file points to correct GitHub repository location
    }
  }
  \subsection{Updates}{
    \itemize{
      \item Added min.cells option for creating Seurat object
      \item Selecting cells in interactive plots will store the cellIds in a variable, for future implementation
      \item Downloading variable gene plot also provides a .txt file with a list of the variable genes found
      \item PC heatmap and Viz plot are on one page, and downloading them also provides a table of genes that make up each PC
      \item 3D plots added for pca and tsne
      \item Point transparency can be changed with the slider input
      \item DEG violin plot is based on a drop-down menu
      \item Tables for DEG list and PC genes improved with the DT package
    }
  }
}

\section{Version 0.0.1 (2017-11-24)}{
  \subsection{Updates}{
    \itemize{
      \item Updated UI to hide PCA plot area until PCA was done
      \item Added FindClusters functionality to PCA tab
      \item PCA and tSNE co-ordinates will be saved to .txt files together with the PDF plots
      \item Cluster IDs also saved to .txt file when PCA plot is saved
    }
  }
}

\section{Version 0.0.0.9000 (2017-09-27)}{
  \subsection{NEW APP}{
    \itemize{
      \item Shiny App UI implementation of Seurat sc Analysis
      \item Input TPM table generated from fastq alignments to analyse single cell data
      \item Interactive plots with Plotly
      \item Save output plots and tables and pdf and csv respectively, as well as saving seurat object for future analysis
    }
  }
}



