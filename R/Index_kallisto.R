#' Build kallisto index
#' 
#' This function creates a kallisto index, for later use to do read quantification.
#' 
#' @param index.name Specify a filename for the output, filetype .idx will be appended to it
#' @param fasta.file Path to the fasta file that will be used to build the index, can be gzipped
#' @param kallisto.path This argument is to be used if kallisto is not found in the \code{PATH} variable
#' @param specify.kmer.size The length of k-mer to be used in the index. Must be an odd number between 3 and 31
#' @param make.unique Logical argument on whether repeated targets should be given unique names, defaults to false
#' @param output.dir Path to directory to save output index file to
#' @param dry.run If true, will return the command to be passed to kallisto
#' 
#' @return If \code{dry.run = TRUE}, prints the command to be passed to kallisto, otherwise creates an index file at target directory
#' 
#' @export
#' 
#' @author Matthew Myint
#'
kallisto_index <- function(index.name, fasta.file, kallisto.path = "kallisto",
                           specify.kmer.size = NULL, make.unique = FALSE,
                           output.dir = getwd(), dry.run = FALSE){
  ## check that fasta files exist
  if(is.character(fasta.file)){
    if(!file.exists(fasta.file)){
      stop(paste0("File ", fasta.file, " doesn't exist!"))
    }
    # $file.sep or $path.sep
    output_name <- paste0(output.dir, .Platform$file.sep, index.name, ".idx")
    kmer_arg <- "-k 31"
    if(!is.null(specify.kmer.size)){
      if(is.numeric(specify.kmer.size)){
        if(specify.kmer.size >= 3 && specify.kmer.size <= 31){
          if(specify.kmer.size %% 2 != 0){
            kmer_arg <- paste("-k", specify.kmer.size)
          }else{
            stop("kmer size must be an odd integer")
          }
        }else{
          stop("k-mer size must be an odd integer from 3-31")
        }
      }else{
        stop("k-mer selection must be an integer")
      }
    }
    unique <- NULL
    if(make.unique == TRUE){
      unique <- "--make-unique"
    }
    kallisto_arguments <- paste("index", "-i", normalizePath(output_name), kmer_arg, unique, normalizePath(fasta.file, mustWork = NA))
    
    if(dry.run){
      return(paste(kallisto.path, kallisto_arguments))
    }else{
      out <- tryCatch(ex <- system2(kallisto.path, args = kallisto_arguments, stdout = TRUE, stderr = TRUE),
                      warning = function(w){w}, error = function(e){e})
      return(out)
    }
  }else{
    stop("Please provide character vector containing the full path to the fasta file!")
  }
}
