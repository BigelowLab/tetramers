#' Extract a name from a filename
#'
#' @export
#' @param filename character, one or more filenames
#' @return one or more names
input_name <- function(filename = c("foo.fasta", "bar.fasta.gz" )){
  
  
  input_name_one <- function(f){
    ix <- gregexpr(".", f, fixed = TRUE)
    if (ix[[1]][1] < 0) {
      name <- f
    } else {
      n <- length(ix[[1]])
      name <- substring(f, 1, ix[[1]][n]-1)
    }       
    name
  }
 
  sapply(basename(filename), input_name_one)
}

#' Retrieve the normalized 'explanation' of variance for each PC as eigenvalues
#'
#' Returns the 'explanations' of variance as derived from the sdev component of the
#'  PC (principal components)
#'
#' @export
#' @param X the Tetramers R6 class object
#' @param n the number of elements to return, by default the number of PCs
#' @param asVariance if TRUE then use the variance (which is the eigenvalue) instead of the
#     the sdev
#' @return a vector of eigenvalues
getExplanation <- function(X, n = length(X$pick), asVariance=TRUE ){
   p <- X$PC$sdev
   if (asVariance) p <- p^2
   p <- (p/sum(p))[1:n]
   names(p) <- names(X$pick)
   p
}