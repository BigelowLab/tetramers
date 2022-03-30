#' Selects the outliers for PCs 1-nPC using a tidy-approach
#'
#' It's a select-and-remove process.  First order by PCx. Select the outliers and
#' remove contigs with matching names.  Then select a second set of outliers.  So we now have
#' 2-low (negative) and 2-high (positive) outliers along axis PCx.  With this winnowed
#' data set order by PCy. Select the outliers and remove contigs with matching names.  Then 
#' select a second set of outliers along this axis.  Now we have the same kind of double set.
#'
#' @export
#' @param X a Tetramers R6 class reference
#' @return the input X
select_outliers <- function(X){

  # given a table, return the row indoices for first and last rows
  outindex <- function(x){
    c(1, nrow(x))
  }
  
  # Given pairs of PCs, select the extreme picks
  # @param x tibble with [cname, wname, PCn, PCm] where PCn and PCn are PC values for say, PC1 and PC2 (or PC3 and PC4)
  # @param pick the number of outliers to pick from each PC vector
  # @return tibble PC, outlier, cname, wname
  select_pc <- function(x, pick = 2){
  
    orignames <- colnames(x)
    colnames(x) <- c("cname", "wname", "PCA", "PCB")
    y <- x %>%
      dplyr::arrange(.data$PCA)

    # PCx extremes 1
    ix <- outindex(y)
    lohiA1 <- y %>%
      dplyr::slice(ix) 
  
    # remove matching contigs  
    y <- y %>%
      dplyr::filter(!(.data$cname %in% lohiA1$cname))
  
    # PCx extremes 2
    ix <- outindex(y)
    lohiA2 <- y %>%
      dplyr::slice(ix) 
  
    # bind and reorder to lo, lo, hi, hi
    lohiA <- dplyr::bind_rows(lohiA1, lohiA2) %>%
      dplyr::slice(c(1,3,4,2))
    
    # remove matching contigs and sort by PCy
    y <- y %>%
      dplyr::filter(!(.data$cname %in% lohiA$cname)) %>%
      dplyr::arrange(.data$PCB)
  
    # PCy extremes 1
    ix <- outindex(y)
    lohiB1 <- y %>%
      dplyr::slice(ix)
     
    # remove matching contigs  
    y <- y %>%
      dplyr::filter(!(.data$cname %in% lohiB1$cname))
  
    # PCx extremes 2
    ix <- outindex(y)
    lohiB2 <- y %>%
      dplyr::slice(ix) 
  
      # bind and reorder to lo, lo, hi, hi
    lohiB<- dplyr::bind_rows(lohiB1, lohiB2) %>%
        dplyr::slice(c(1,3,4,2))

    lohiA <- lohiA %>%
      dplyr::select(.data$cname, .data$wname) %>%
      dplyr::mutate(
            PC = orignames[3],
            pick = c("lo", "lo", "hi", "hi"))
    lohiB <- lohiB %>%
      dplyr::select(.data$cname, .data$wname) %>%
      dplyr::mutate(
            PC = orignames[4],
            pick = c("lo", "lo", "hi", "hi"))
    
    dplyr::bind_rows(lohiA, lohiB)
  } #select_pc
    

  pick <- X$pick
  nPC <- length(pick)
  M <- X$PC$x[,1:nPC, drop = FALSE]
  wname <- rownames(M)
  cname <- as.vector(decomposeContigNames(wname)[,1])
  ucname <- unique(cname)
  nucname <- length(ucname)
  if (nucname < nPC){
    if (is_odd(nucname)){
      nPC <- nucname - 1
    } else {
      nPC <- nucname
    }
    if (nPC <= 1) stop("unable to proceed with just one contig")
  }
  
  
  x <- dplyr::tibble(cname = cname, wname = wname) %>%
    dplyr::bind_cols(dplyr::as_tibble(M))
  index <- seq(from = 1, to = nPC, by = 2)
  lapply(index,
    function(i){
      pcs <- paste0("PC", c(i, i+1))
      y <- select_pc(x %>% dplyr::select(cname, wname, pcs), 
                     pick = pick[[pcs[1]]])
    }) %>%
    dplyr::bind_rows()
}