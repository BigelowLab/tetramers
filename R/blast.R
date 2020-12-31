#' Run blast and extract the results
#'
#' @export
#' @param X a Tetramers R6 class object
#' @param verbose logical, if TRUE output messages for debugging purposes
#' @return the updated input object, blast results stored in the 'blast_data' field
#'  and blast info stored in 'blast_info' field
run_blast_app <- function(X, verbose = FALSE){
  verbose <- verbose[1]
  if (!inherits(X, "Tetramers")) stop("input must be Tetramers R6 class object")
  
  blastOpts <- strsplit(X$blast_options, " ")[[1]]
  blastApp <- Sys.which(blastOpts[1])
  if (nchar(blastApp) == 0){
      cat("blast is not accessible\n")
      return(invisible(X))
  }
  query_file = file.path(X$out_dir, 
      paste(X$name,"outliers.fasta",sep="-"))
  xml_file = file.path(X$out_dir, 
      paste(X$name,"outliers.xml",sep="-"))

  CMD <- paste(X$blast_options, "-outfmt 5", 
      "-query", query_file, "-out", xml_file)
  if (verbose) cat("running...", CMD, "\n")
  ok <- system(CMD)
  if (ok != 0){
       cat("blast command failed\n")
       cat(CMD, "\n")
       return(invisible(X))
  }
     
  info = try(blastxml::blastxml_info(xml_file))
  if (!inherits(info, 'try-error')) X$blast_info <- info
  
  x <- try(blastxml::blastxml_dump(xml_file, form = 'tibble'))
  if (!inherits(x, 'try-error')) X$blast_data <- x
  
  invisible(X)
}

#' Trim a blast query name to a reasonable length
#' 
#' @param x blast tibble
#' @param len numeric, the maximum string length allowed
#' @return a character vector
trim_QueryDef <- function(x, len = 30) {

   if (is.list(x)) {
      txt <- lapply(x, trim_QueryDef, len = len)
   } else {
      txt <- sprintf(paste0("%0.", len,"s"), x[["Hit_def"]])     
   }
  
  txt
}


#' Trim a blast hit description string to a particular length with it's associated
#' hit score and accession name.
#' 
#' @param x blast tibble
#' @param len numeric, the maximum output string length
#' @param hsp_bit_score_min numeric, hits below this are flagged with noHitText
#' @param noHitText the value to return for failed hits
#' @return the tibble with an added column 'hit_tetxt' 
blastHitText <- function(x, len = 30, 
   hsp_bit_score_min = 75,
   noHitText = "no significant hits"){

   txt <- paste(
     sprintf("%0.10s bit score = %0.0f", 
             x[["Hit_accession"]],  
             x[["Hsp_bit-score"]]),
     sprintf(paste0("%0.", len,"s"), 
             gsub("[\t\n]", " ", x[["Hit_def"]])), 
     sep = "\n")

   txt[(x[['Hit_len']] <= 0) | (x[['Hsp_bit_score']] <= hsp_bit_score_min)] <- noHitText  
   names(txt) <- x[['Iteration_query-def']]
    
   txt
}
