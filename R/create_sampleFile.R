#' Create Sample Files
#' 
#' This function creates a file telling the package which files to obtain sequence data from.
#' Depending on the user's input, the output will be catered to single or paired read experiments
#' 
#' @param fileName1 For single-end reads, a vector of all fastq files. For paired end, a vector of all fastq that are first in pair
#' @param fileName2 Vector of filepaths for all fastq for reads that are second in their paired reads. If empty, will treat files in \code{fileName1} as single-end
#' @param sampleNames Vector of names for samples
#' @param write.to.file Logical indicating if the table should be written to a file, which is used for the quantification function
#' @param output.name User-specified filename for the output file
#' 
#' @return A table with either 2 (single-end) or 3 (paired-end) columns
#' 
#' @note Ensure that the number of reads matches up with the number of sample names given. Also, due to
#' constraints of \code{runKallisto()}, the sample file created can only be saved to where the fastq files are.
#' It is also necessary to keep all fastq files in the same directory
#' 
#' @export
#' 
#' @author Matthew Myint
#' 
createSampFile <- function(fileName1, fileName2, sampleNames, write.to.file = FALSE, output.name = "sampleFiles"){
  ## Argument checks
  if(missing(fileName1)){
    stop("No files provided, please provide a vector of filepaths!")
  }
  ## for single end reads
  if(missing(fileName2)){
    if(missing(sampleNames)){
      warning("No sample names provided, using filenames as sample names!")
      sampleNames <- sub("\\..+", "", basename(fileName1))
    }
    sampFile <- list(SampleName = sort(sampleNames), FileName = sort(fileName1))
    sampFile <- as.data.frame(sampFile)
  }else{
    ## Paired reads
    ## check that reads are actually paired
    if(!length(fileName1) == length(fileName2)){
      stop("Length of filenames 1 and 2 do not match. Please check your input files!")
    }
    if(missing(sampleNames)){
      warning("No sampleNames provided, using filenames as sample names!")
      sampleNames <- sub("(_[0-9]\\.).+", "", basename(fileName1))
    }
    check <- list(fileName1, fileName2)
    check <- lapply(check, function(x) sub("(_[0-9]\\.).+", "", basename(x)))
    if(!length(setdiff(check[[1]],check[[2]])) == 0){
      stop(paste0("The following files may not be paired. Please amend filenames as needed.", "\n",
                  "From FileName1: ", paste(setdiff(fileName1, fileName2), collapse = " "), "\n",
                  "From FileName2: ", paste(setdiff(fileName2, fileName1), collapse = " ")))
    }
    sampFile <- list(SampleName = sort(sampleNames),
                     FileName1 = basename(sort(fileName1)),
                     FileName2 = basename(sort(fileName2)))
    sampFile <- as.data.frame(sampFile)
  }
  ## due to constraints of runKallisto(), force sampleFile to be written to same directory as fastq files
  fastq.Dir <- dirname(fileName1[1])
  
  ## write to text
  if(write.to.file){
    write.table(sampFile, file = paste0(fastq.Dir, "/", output.name, ".txt"), sep = "\t", row.names = FALSE)
  }
  return(sampFile)
} 
