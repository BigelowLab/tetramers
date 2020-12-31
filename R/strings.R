#' Tidy an input FASTA string
#'
#' @export
#' @param x the input contig sequence OR a simple string
#' @param dict character, the dictionary of allowed characters
#' @param replacement character, the replacement value for disallowed characters
#' @return either a contig squence or a character string - matches the input
clean_sequence <- function(x, 
  dict = Biostrings::DNA_ALPHABET, 
  replacement = "N"){

    oldx <- paste0(tolower(dict), collapse = "")
    newx <- paste0(dict, collapse = "")
    x <- Biostrings::chartr(oldx, newx, x)

    uletters <- Biostrings::uniqueLetters(x)
    ok <- !(uletters %in% dict)

    if (any(ok)) {
        oldx <- uletters[ok]
        newx <- rep(replacement, length(oldx))
        x <- Biostrings::chartr(oldx, newx, x)
    }
    x
}

#' Parse the name of a contig
#'
#' Given the contig-window names "blahblah_1-5000" return a 2 element
#' character vector of c("blahblah", "1-5000")
#'
#' @export
#' @param x character vector of contig names
#' @param sep character separator around which the component names are split
#' @return a nrowx2 character array
decomposeContigNames <- function(x, sep = "_"){
   ix <- gregexpr("_", x)
   last_ <- function(x) return(x[length(x)])
   lastIx <- sapply(ix, last_)
   split_ <- function(x, i) substring(x, c(1,i+1), c(i-1,nchar(x)))
   t(mapply(split_, x, lastIx))
}


#' Retrieve the window specification from a window name.
#'
#' @export
#' @param X a TetramerPCA object
#' @param windowName the full name of the window ("contigName_start-end")
#' @param delim the delimiter between the start and end of the window
#' @return the sequence within the window
getBlastWindow <- function(X, windowName = "unknown_1-5000", delim = "-"){
   s <- decomposeContigNames(windowName)
   ss <- as.numeric(unlist(strsplit(s[2], delim)))
   as.character(X$rawdata[[s[1]]][ss[1]:ss[2]])
}

#' Create a  vector of the default combinations of A, C, G and T
#'
#' @export
#' @param what character, the name of the combinations, "tetramers" means A,C,G, and T
#'    while "tetramers-all" adds "N" to the tetramers list.
#' @param sorted logical, sort into lexical order (AAAA, AAAC, ..., TTTT)
#' @param as.DNAStringSet logical, if TRUE then return a DNAStringSet otherwise a vector
#' @param lowercase logical, if TRUE then use lower case letters
#' @return a character vector or DNStringSet
makeCombinations <- function(what = c("tetramers", "tetramers-all")[1],
   sorted = TRUE, as.DNAStringSet = FALSE, lowercase = FALSE){

   if (lowercase) {
      p <- switch(what,
         "tetramers" = c("a", "c", "g", "t"),
         "tetramers-all" = c("a", "c", "g", "t", "n"))
   } else {
      p <- switch(what,
         "tetramers" = c("A", "C", "G", "T"),
         "tetramers-all" = c("A", "C", "G", "T", "N"))
   }
   p <- as.matrix(expand.grid(p,p,p,p, stringsAsFactors = FALSE))
   p <- apply(p, 1, "paste", collapse = "")
   if (sorted) p <- sort(p)
   if (as.DNAStringSet) {
      p <- Biostrings::DNAStringSet(p)
   }
   return(p)
}

