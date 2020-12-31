#' Tabulate the tetramers within a window advanced by step
#'
#' This generally is called with a DNAStringSet, and then is called recursively
#' with DNAString.
#'
#' @export
#' @param x a DNAStringSet (or DNAString)
#' @param tWindow numeric, the width of the sliding window in base units
#' @param tStep numeric, the number of bases to advance between successive window locations
#' @param tWidth numeric, the width of the genome unit, for tetramer analysis it's 4
#' @return for each DNAString a list of \code{status} (0/1 where 0 = ok) and \code{x} a
#'    tabulation matrix.  In the matrix each row corresponds to a window location along
#'    the length of the sequence.
tabulate_windows <- function(x,
   tWindow = 1600, tStep = 200, tWidth = 4){
  
  DNALETTERS <- TETRA[["LETTERS"]]
  UTETRAMERS <- TETRA[["UTETRAMERS"]]
  nutet <- length(UTETRAMERS)
  tetnm <- names(UTETRAMERS)
  N <- length(x)
  starts <- seq.int(from = 0, to = tWindow - tWidth)
  nstarts <- length(starts)

  if (inherits(x, "XString")){

    DNALETTERS <- TETRA[["LETTERS"]]
    UTETRAMERS <- TETRA[["UTETRAMERS"]]
    nutet <- length(UTETRAMERS)
    tetnm <- names(UTETRAMERS)
    N <- length(x)
    starts <- seq.int(from = 0, to = tWindow - tWidth)
    nstarts <- length(starts)

    if ( N >= tWindow ) {
       nstep <- ( ( N-tWindow ) %/% tStep ) + 1
       M <- matrix(0, nrow = nstep, ncol = nutet)
       colnames(M) <- tetnm
       nm <- vector(mode = "character", length = nstep)

       for (istep in seq_len(nstep)){
          #compute where we are along the contig sequence
          s0 <- (istep-1) * tStep + 1
          # compute ends from reusable starts
          s1 <- s0 + starts
          #make the name
          nm[istep] <- sprintf('%i-%i', s1[1], s1[nstarts] + tWidth - 1)
          # split the string into tetramers
          s <- as.character(Biostrings::DNAStringSet(x, start = s1, width = tWidth))
          # filter for only combinations of AGCT, stuff those into the matrix
          ix <- match(s, tetnm, nomatch = 0)
          M[istep, ] <- tabulate(ix, nbins = nutet)

       } # istep-loop

       rownames(M) <- nm
       tab <- list(status = 0, x = M)

    } else {  # now what to do with short ones .... N < tWindow
       # generate a table of the bases present
       tab <- Biostrings::letterFrequency(x, names(DNALETTERS))
       tab <-  list(status = 1, x = tab)

    }

  } else if (inherits(x, "XStringSet")) {

     tab <- lapply(x, tabulate_windows, tWindow = tWindow,
                   tStep = tStep, tWidth = tWidth)
  }

  return(tab)
}


#' Reduce the tetra census table to include just the 136 complements (for tetramers)
#' by pooling reverse complements.
#'
#'
#' @export
#' @param x a list of matrices or a single matrix (nWindow x 256)
#' @param comp a named character vector of reverse complements
#' @param keep a character vector of name tetramers to retain
#' @return a list or matrix with just 136 columns
reduce_tetra <- function(x, 
  comp = TETRA[["RCTETRAMERS"]],
  keep = TETRA[["TETRAMERS"]]){
  if (is.matrix(x)){
     nm <- names(comp)
     x[,nm] <- x[,nm, drop = FALSE] + x[,comp, drop = FALSE]
     x <- x[,keep, drop = FALSE]
  } else {
     x <- lapply(x, reduce_tetra, comp = comp, keep = keep)
  }
  return(x)
}

#' Normalize the tetra sequences
#'
#' This function normalizes the tetramer counts by window locations
#' You can visualize each window location by one row of the census matrix
#' Normalization occurs across each row, with the max count scoring the highest
#' value toward 1.0 and the least represented as toward 0.0.
#' If there are no hits of the any TETRAMERS (all Ns in the original)
#' then the normalized score for each tetramer in the row is 0.
#'
#' @export
#' @param x a list of matrices or a single matrix
#' @return a list or matrix of the same form as the input
normalize_tetra <- function(x){
   if (is.matrix(x)) {
      s <- rowSums(x, na.rm = TRUE)
      ix <- (s <= 0)
      x <- sweep(x, 1, s, "/")
      x[ix,] <- 0
   } else {
      x <- lapply(x, normalize_tetra)
   }
   return(x)
}

#' Bind a list of tetra matrices to one matrix (n x 136)
#'
#' @export
#' @param x a list of tetramer matrices
#' @return a matrix with uniquely named rows
bind_tetra <- function(x){
   nm <- names(x)
   for (n in nm)  rownames(x[[n]]) <- paste0(n, "_", rownames(x[[n]]))
   do.call(rbind, x)
}

#' Perform Principal Component Analysis on the normalized, reduced, tetramer matrix
#'
#' @export
#' @param x numeric, the (n x 136) matrix of normalized counts
#' @param scale. see \code{\link{prcomp}}
#' @param ... further arguments for \code{\link{prcomp}}
pca_tetra  <- function(x, scale. = TRUE, ...){
      xcols <- apply(x, 2, sum)
      x <- x[,xcols > 0]
      PC <- try(prcomp(x, scale. = scale., ...))
      if (inherits(PC, "try-error")){
         cat(str(x), "\n")
         return(list())
      }
      class(PC) <- list("list", "prcomp")
      invisible(PC)
}
