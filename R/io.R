#' Read a FASTA file and possibly tidy it
#'
#' @export
#' @param X Tetramers-class object
#' @param clean logical, if TRUE then tidy the input
#' @return X Tetramers-class object
read_fastafile <- function(X, clean = TRUE){

    if (missing(X)) stop("Input Tetramers-class object is required")

    x <- try(Biostrings::readBStringSet(X$filename[1]))
    if (inherits(x, 'try-error')){
        cat('error reading fasta file\n')
        return(X)
    }

    if (clean) {
         x <- clean_sequence(x)
    } else {
        oldx <- paste0(tolower(Biostrings::DNA_ALPHABET), collapse = "")
        newx <- paste0(Biostrings::DNA_ALPHABET, collapse = "")
        x <- Biostrings::chartr(oldx, newx, x)
    }

    X$rawdata <- Biostrings::DNAStringSet(x)
    invisible(X)
}


#' Read a blast xml output file
#' 
#' @export
#' @param filename character, the name of the xml file to read
#' @return a lit with two elements: info and data either of which may be NULL
read_blast_output <- function(filename){
  stopifnot(file.exists(filename[1]))
  info <- try(blastxml::blastxml_info(filename[1]))
  if (inherits(info, 'try-error')) info <- NULL
  
  x <- try(blastxml::blastxml_dump(filename[1], form = 'tibble'))
  if (inherits(x, 'try-error')) x <- NULL
  
  invisible(list(info=info, data = x))
}


#' Write all of the known results
#'
#' @export
#' @param X Tetramers-class object
#' @param path the output path - created if it doesn't exist
#' @param what character - a list of items to write
#' \itemize{
#'  \item{input the input FASTA file if it is available and doesn't already exist in path}
#'  \item{info a dump of input parameters to a yaml}
#'  \item{count tetramer counts}
#'  \item{pca pca values}
#'  \item{fail listing of failed tetramer windows}
#'  \item{loading pca loading values}
#'  \item{outliers table of outliers}
#'}
#' @param verbose logical, if TRUE then output messages
write_tetra <- function(X,
  path = X$out_dir,
  what = c("input", "info", "count", "pca", "fail", "loading", "outliers"),
  verbose = FALSE){
  
  what <- tolower(what)
  verbose <- verbose[1]
  path <- path[1]
  if (!dir.exists(path)){
    if (verbose) cat("creating output directory\n")
    ok <- dir.create(path, recursive = TRUE)
    if (!ok) stop(paste("Unable to create the output path:", path))
  }
  
  if ("input" %in% what){
    
    filename <- file.path(path, basename(X$filename))
    if (file.exists(X$filename) & !file.exists(filename)){
      if (verbose) cat("copying input FASTA to output\n")
      ok <- file.copy(X$filename, path, overwrite = TRUE)
    } else {
      if (verbose) cat("input FASTA already exists in output directory\n")
    }
  }
  
  if ("info" %in% what){
    if (verbose) cat("saving metadata in yaml\n")
    filename <- file.path(path, paste(X$name,"info.yaml",sep="-"))
    info <- yaml::as.yaml(list(
      input = X$filename,
      name = X$name,
      out_dir = X$out_dir, 
      pick = as.list(X$pick),
      parameters = as.list(X$parameters),
      blast_options = X$blast_options,
      R_version = R.version.string))
      
    cat(info, sep = "", file = filename)
  }
  
  if ("count" %in% what){
    if (verbose) cat("saving tabulation\n")
    filename <- file.path(path, paste(X$name,"tetramer-counts.csv.gz",sep="-"))
    readr::write_csv(X$get_tabulation(form = 'table'), filename)
  }
  
  if ("pca" %in% what){
    if (verbose) cat("saving principal components\n")
    filename <- file.path(path, paste(X$name,"tetramer-PC.csv.gz",sep="-"))
    readr::write_csv(X$get_principal_components(form = 'table'), filename)
  }
  
  if ("fail" %in% what){
    if (verbose) cat("saving fails (gzipped if any)\n")
    filename <- file.path(path, paste(X$name,"tetramer-fail.csv.gz",sep="-"))
    if (length(X$fail) > 0){
       z <- do.call(rbind,lapply(X$fail, "[[", "x"))
       status <- sapply(X$fail, "[[", "status")
       contigLength <- apply(z,1,sum)
       nz <- z/contigLength
       colnames(nz) <- paste("p", colnames(z), sep = "")
       pGC <- nz[,"pC"] +  nz[,"pG"]
       dplyr::tibble(status = TETRA[["FAILSTATUS"]][status], 
                          contigLength = contigLength, 
                          pGC = signif(pGC,digits = 3)) %>%
        readr::write_csv(filename)
    } else {
      cat("No failed contigs\n", 
        file = file.path(path, paste(X$name,"tetramer-fail.csv",sep="-")))
    }
  }
  
  if ("loading" %in% what){
    if (verbose) cat("saving loadings\n")
    filename <- file.path(path, paste(X$name,"tetramer-loading.csv.gz",sep="-"))
    readr::write_csv(X$get_loadings(form = 'table'), filename)
  }
  
  if ("outliers" %in% what){
    # table of info
    if (verbose) cat("saving outlier info table\n")
    filename <- file.path(path, paste(X$name,"outliers.csv",sep="-"))
    readr::write_csv(X$outliers, filename)

    # now we need a FASTA for the outlier window sequences
    # a little function to write one contig
    # @param x one row tibble from get_outlier_sequences()
    # @param y ignored, for futire use with group_walk()
    # @param handle file handle
    # @param width numeric, the number of columns before wrapping
    # @return nothing
    write_contig <- function(x, y,  handle = NULL, width = 80){
      cat(sprintf(">%s", x$wname),"\n", sep = "", file = handle)
      nx <- nchar(x$seq)
      ix <- seq(from = 1, to = nx, by = width)
      iy <- seq(from = width, to = nx, by = width)
      if (iy[length(iy)] < nx) iy <- c(iy, nx)
      xs <- Biostrings::DNAStringSet(x$seq, start = ix, end = iy)
      cat(as.character(xs), sep = "\n", file = ff)
    }
    if (verbose) cat("saving outlier fasta file\n")
    filename <- file.path(path, paste(X$name,"outliers.fasta",sep="-"))
    ff <- file(filename, open = "wt")
    x <- X$get_outlier_sequences()
    for (i in seq_len(nrow(x))){
      write_contig(x %>% dplyr::slice(i), handle = ff)
    }
    close(ff)
    
  }
  
  invisible(NULL)
}
