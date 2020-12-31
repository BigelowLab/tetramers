#' Compute pretty [xy] ranges for TetramerPCA data
#'
#' @export
#' @param X A TetramerPCA class object
#' @param n the number of ranges to return, by default the length of PCs we pick
#' @return a 2 column by n row matrix of min and max ranges
getTetramerRange <- function(X, n = length(X$pick)){
      #            min       max
      # PC1 -20.407072 12.167675
      # PC2  -9.930063  8.536594
      # PC3  -8.646382 28.800613
      # PC4  -8.825642 10.554186
      # PC5  -8.765708  7.011479
      # PC6 -55.177973  5.763522
      # PC7  -6.080349 15.275533
      # PC8 -46.890103  9.361013
   xmax <- apply(X$PC$x, 2, max)
   xmin <- apply(X$PC$x, 2, min)
   cbind(min = xmin[1:n], max = xmax[1:n])
}


#' trims input text vector to length and adds etc if needed in such a way that
#' \code{nchar(paste(newX, etc)) < length}
#'
#' @export
#' @param x character vector of one or more strings
#' @param length numeric, maximum string length allowed
#' @param etc character, the flag to indicate removed characters (a suffix)
#' @param strip.white logical, if TRUE strip newlines and tabs
#' @return input chacater vector with elements possibly abbreviated
trimLegendText <- function(x, length = 35,
        etc = "...", strip.white = TRUE){
   if (strip.white) x <- gsub("[\t\n]", " ", x)
   len <- sapply(x, nchar)
   fixMe <- len > length
   x[fixMe] <- paste(substring(x[fixMe], 1, length-nchar(etc) - 1), etc)
   x
}

#' Returns a numeric vector suitable of pch values for plot symbols
#' 
#' @export
#' @param x a numeric vector - we only care about its length
#' @param what the symbol identifier
#' @return the value for what repeated
pchfun <- function(x = 1:10, what=24) rep(what, length(x))

#' Returns a numeric vector suitable of cex values for plot symbols
#'
#' @export
#' @param x a numeric vector - we only care about its length
#' @param big numeric, cex values for big symbols
#' @param small numeric, cex values for small symbols
#' @param what character
#' \itemize{
#'   \item{'simple' are small size}
#'   \item{'linear' points are scales from big to small}
#'   \item{'bigfirst' first point is big all others small}
#' }
#' @return numeric vector suitable for cex
cexfun <- function(x=1:10, big = 1.5, small = 0.5,
    what = c('simple', 'bigfirst', 'linear')[1]){
    switch(what[1],
        'linear' = seq(from = big, to = small, length = length(x)),
        'bigfirst' = c(big, rep(small, length(x)-1)),
        rep(small, length(x)))
}


#' Returns the even-numbered color palette associated with the count provided.
#'
#' @param x - the number of color steps, if odd then n+1 is returned,
#'  if n > maxN (see \code{RColorBrewer}) the maxN colors are returned with a
#'  warning otherwise n colors are returned
#' @param name the name of the ColorBrewer palette to return
#' @return the ColorBrewer palette specified
get_pal <- function(x, name = "PuOr"){
   x = x + ifelse((x %% 2) == 1, 1, 0)
   RColorBrewer::brewer.pal(x, name)
}

#' Plot one PCx vs PCy pairing
#'
#' @export
#' @param X a TetramerPCA class object
#' @param PC two element vector of PC names
#' @param MAXSTRING numeric, the maximum displayed length of a contig name
#' @param noHitText the value to return for failed hits
#' @param addBlast, logical
plot_PCvPC <- function(X,
    PC = c("PC1", "PC2"),
    MAXSTRING = 35,
    noHitText = 'no significant hits',
    addBlast = !is.null(X$blast_data)){

    if (FALSE){
      PC = c("PC1", "PC2")
      MAXSTRING = 35
      noHitText = 'no significant hits'
      addBlast = !is.null(X$blast_data)
    }

    # what a mess of base-graphics startup stuff
    par(omd = c(0.05, 0.95, 0.15, 1))
    lim <- getTetramerRange(X)
    pcex = c(1.5, 0.4)
    px <- PC[[1]]
    py <- PC[[2]]
    xlim <- lim[px,] + c(-10,10)
    ylim <- lim[py,] + c(-3,3)
    ex <- matrix(getExplanation(X)[PC]*100 , byrow = TRUE, ncol = 2)
    plabel <- sprintf("%s, %0.1f%% variation explained", PC, ex)
    names(plabel) <- PC
    def <- X$parameters
    tetraWindow <- def[["WINDOW"]]
    tetraStep <- def[["STEP"]]
    scoreName <- "Hsp_bit-score"
    
    # here we grab the PC values but we toss all but names and PCx and PCy
    PCV <- X$get_principal_components() %>%
      dplyr::select(.data$cname, .data$wname, PC)
    
    # here we plot all PCx and PCy points as grey dots
    plot(PCV[[px]], PCV[[py]],
        main = X$name,
        xlim = xlim, 
        ylim = ylim,
        col = "grey75",
        cex = 0.3,
        pch = 16,
        xlab = plabel[1],
        ylab = plabel[2])
    # provode user with recall nudge about windowing parameters
    mtext(sprintf("window = %i step = %i", tetraWindow, tetraStep),
            side = 3, line = 0.8, cex = 0.7)
    
    # Given one row of the OUT table
    # plot the points from the PCV table that match the 
    # contig name (cname)
    # adds a nwindow count to returned table (same as input with one new variable)
    # 
    # as a side effect, if dev.cur() > 1 then points are plotted
    #
    # @param out one row tibble sliced from OUT
    # @param y ignored
    # @param px and py character vectors of the PC columns
    # @return input one-row table with nwindow added
    #   note the 
    plot_window <- function(out, y, px = "PC1", py = "PC2"){
      x <- PCV %>%
        dplyr::filter(.data$cname %in% out$cname)
      if (dev.cur() > 1){
        points(x[[px]], x[[py]],
           pch = out$pch,
           col = out$col,
           bg =  out$col,
           cex = cexfun(x[[px]]),
           typ = 'o')
      }
      out %>%
          dplyr::mutate(nwindow = nrow(x))
    }

    plotch <- c(hi = 24, lo = 25)
    n = X$pick[[px]] * 2
    cols <- c(get_pal(n, "PuOr"), get_pal(n, "PiYG"))

    OUT <- X$outliers %>%
      dplyr::filter(.data$PC %in% c(px, py)) %>%
      dplyr::mutate(pch = plotch[.data$pick],
                    col = cols)
      # # here, using an updated dplyr we might use group_map
      # # but in this old version (pre v1) we have to lapply
      # dplyr::rowwise() %>%
      # dplyr::group_map(plot_window, px = px, py = py, .keep = TRUE) %>%
      # dplyr::bind_rows()    
    OUT <- lapply(seq_len(nrow(OUT)),
      function(i){
        OUT %>%
          dplyr::slice(i) %>%
          plot_window(y = NULL, px = px, py = py)
      }) %>%
      dplyr::bind_rows()
                    
    ###
    # finally the legend
    ###
    # X-legend on left
    out <- OUT %>% 
      dplyr::filter(.data$PC %in% px) %>%
      dplyr::slice(4:1)
    ix <- !is.na(out$cname)
    if (any(ix)){
       lgn <- sprintf(paste0("%0.",MAXSTRING,"s (%0.0d)"),
           trimLegendText(out$cname[ix]), out$nwindow[ix])
       u <- par("usr")
       xL <- u[1]
       yT <- grconvertY(0.05, from ="nfc", to = "user")

       legend(xL, yT, lgn, title = px, title.col = "black",
          pch = out$pch,       #c(phi1,plo1),
          pt.bg = out$col,     #c(chi1,clo1),
          col = out$col,       #c(chi1,clo1),
          text.col = out$col,  #c(chi1,clo1),
          xpd = NA, bty = "n", cex = 0.7, xjust = 0)
    }
    # Y-legend on right
    out <- OUT %>% 
      dplyr::filter(.data$PC %in% py) %>%
      dplyr::slice(4:1)
    ix <- !is.na(out$cname)
    if (any(ix)){
       lgn <- sprintf(paste0("%0.",MAXSTRING,"s (%0.0d)"),
           trimLegendText(out$cname[ix]), out$nwindow[ix])
       xL <- u[2]
       legend(xL, yT, lgn, title = py, title.col = "black",
          pch = out$pch,       #c(phi1,plo1),
          pt.bg = out$col,     #c(chi1,clo1),
          col = out$col,       #c(chi1,clo1),
          text.col = out$col,  #c(chi1,clo1),
          xpd = NA, bty = "n", cex = 0.7, xjust = 1)
    }
      
    
   
   if (addBlast){
     # provide use with info about blast if applicable
     # note date is date graphics run not the date blast was run (which is often the same)
       hsp_bit_score_min = X$parameters[['HSP_BIT_SCORE_MIN']]  #75
       bst = if (hsp_bit_score_min > 0){
           sprintf(" bit score threshold %i,", hsp_bit_score_min )
       } else {
           ""
       }
       mtext(
         sprintf("%s against GenBank %s,%s %s",
                 X$blast_info[['BlastOutput_version']], 
                 X$blast_info[['BlastOutput_db']], 
                 bst,
                 format(Sys.time(), "%Y-%m-%d")), 
         side = 3, line = 0.2, cex = 0.7)
       blast.cex = 0.5
   
       #dx = -0.1    # how far left/right to adjust text alignment

       oxpd <- par("xpd")
       par(xpd = NA)
       
       # first get a vector of the hit_text, but keep in mind that 
       # not every contig window generates a blast result, so we have to match by
       # contig names and window name as needed
       blast_data <- X$blast_data
       blast_data$hit_text <-  unname(blastHitText(blast_data ,
                                      hsp_bit_score_min = hsp_bit_score_min,
                                      noHitText = noHitText))
       
       # split blast_data into window groups, and 
       # then take the first hit for each. I don't know why first
       # but it was done from early on
       # when we get to update to dplyr v1+ we could make all this prettier
       ff <- blast_data[['Iteration_query-def']]
       blast_data <- blast_data %>%
         split(ff) %>%
         lapply(head, 1) %>%
         dplyr::bind_rows() %>%
         dplyr::select(wname = 'Iteration_query-def',
                       hit_text = .data$hit_text) %>%
         dplyr::filter(.data$wname %in% OUT$wname)
       
       # essentially a left join
       OUT$hit_text <- noHitText
       ix <- match(OUT$wname, blast_data$wname)
       OUT$hit_text <- blast_data$hit_text[ix]
       
       # now for each OUT window, find the extreme x and y values
       # and noodle out text position relative to point

       ix <- match(OUT$wname, PCV$wname)
       OUT <- OUT %>%
        dplyr::mutate(x = PCV[[px]][ix],
                      y = PCV[[py]][ix])
        
        # I hate this part adj is [left-right, up-down]
        # "Values of 0, 0.5, and 1 specify that (x, y) should align with the left/bottom, 
        # middle and right/top of the text, respectively."
        # but we have left most, right most, bottom most and top most
        dx <- c(lo = 1, hi = 0) 
        dy <- c(lo = 1, hi = 0)
        pos <- c(lo = 2, hi = 4)
        OUT <- split(OUT, OUT$PC) %>%
        lapply(function(x){
          if (x$PC[1] == px) {
            x$pos <- pos[x$pick]
          } else {
            x$pos <- 4
          }
          x
          }) %>% 
        dplyr::bind_rows()
      
        text(OUT$x, OUT$y, OUT$hit_text,
          col = "blue",
          cex = blast.cex,
          pos = OUT$pos)
   } # addBlast?
   
    
}

#' Plot the Tetramers results with our without blast annonations
#'
#' @export
#' @param X a Tetramers class object
#' @param filename character, the name of the output file
#' @param dim a 2 element vector of output width and height in inches
#' @param MAXSTRING numeric, the maximum displayed length of a contig name
#' @param addBlast logical for suppressing blast annotations
#' @param ... further arguments for the \code{pdf} device
#' @return named logical indicating that the ouptut file exists
plot_tetramers <- function(X,
    filename = file.path(X$out_dir,
        paste0(paste(X$name ,X$parameters[["WINDOW"]],
            X$parameters[["STEP"]], "PC",sep="-"),".pdf")),
    dim = c(width = 8, height = 8),
    MAXSTRING = 35,
    addBlast = !is.null(X$blast_data),
    ...){

    if (!inherits(X,"Tetramers"))
            stop("Input must be Tetramers R6 class object")

    pdf(filename, width = dim[[1]], height = dim[[2]], ...)
      Pairs = matrix(names(X$pick), byrow = TRUE, ncol = 2)
      for (i in 1:nrow(Pairs)){
        plot_PCvPC(X, 
                   PC = Pairs[i,], 
                   addBlast = addBlast,
                   MAXSTRING = MAXSTRING)
        } # i-loop through pairs
    dev.off()
      
    invisible(file.exists(filename))
}