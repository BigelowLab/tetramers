#' A base container class for tetramer analysis
#'
#' @description R6 base class to hold all of the goodies
#' @export
Tetramers <- R6::R6Class("Tetramers",
  public = list(
    #' @field filename character, the fully qualified FASTA filename
    filename = NULL,
    #' @field name character, parsed from the filename by default but settable
    name = NULL,
    #' @field out_dir character, the fully qualified output path, possibly created
    out_dir = NULL,
    #' @field parameters named numeric vector, composed of 
    #'  \itemize{
    #'    \item{window window size in bases, default = 1600}
    #'    \item{window step the distance in bases to advance the window,
    #'      default = 200}
    #'    \item{width sort of dumb, but tetramers are in 4-bases, default = 4}
    #'    \item{hsp_bit_score_min hits below this are flagged with noHitText,
    #'      default = 75}
    #'  }
    parameters = NULL,
    #' @field pick numeric, named vector of number out outliers to pick per PC - generally 2.
    pick = NULL,
    #' @field blastOpts list, a named list of blast options
    blast_options = NULL,
    #' @field blast any, blast results
    blast_data = NULL,
    #' @field blastinfo any, info on the blast application
    blast_info = NULL,
    #' @field rawdata any, raw counts
    rawdata = NULL,
    #' @field data matrix, 
    data = NULL,
    #' @field fail list, list of failed contigs (generally too short)
    fail = NULL,
    #' @field PC list, prcomp results of which we are interested in "x"
    PC = NULL,
    #' @field outliers tibble, selected outliers per PC
    outliers = NULL,

    #' @description Create a new Tetramers container object
    #' @param x character, filename for FASTA format sequences
    #' @param name character, by default extracted from filename
    #' @param output_dir the name of the output directory, if not provided then it is automatically a
    #'    subdirectory of the input file's directory
    #' @param parameters numeric, named vector of sleection parameters
    #'  \itemize{
    #'    \item{window window size in bases, default = 1600}
    #'    \item{window step the distance in bases to advance the window,
    #'      default = 200}
    #'    \item{width sort of dumb, but tetramers are in 4-bases, default = 4}
    #'    \item{hsp_bit_score_min hits below this are flagged with noHitText,
    #'      default = 75}
    #'  }
    #' @param pick numeric, named numeric vector of PC id with number to pick.
    #'   By default 2 extremes  (2 pos, 2 neg) for PC1-8.
    #' @param select_method character, by default "tidy". Ignored.
    #' @param blast_options character, the blast options
    initialize = function(x,
      name = input_name(x[1]),
      output_dir = file.path(dirname(x[1]), name), 
      parameters = c(WINDOW = 1600, STEP = 200, WIDTH = 4, HSP_BIT_SCORE_MIN = 75),
      pick = c(PC1 = 2, PC2 = 2, PC3 = 2, PC4 = 2, PC5 = 2, PC6 = 2, PC7 = 2, PC8 = 2),   
      blast_options = "blastn -db nt -num_threads 4 -num_alignments 10 -evalue 10"){
        
     
      if (!file.exists(x[1])) stop("input filename must be provided")
      self$filename       <- x[1]    
      self$name           <- name[1]
      self$out_dir        <- output_dir
      self$parameters     <- parameters
      self$pick           <- pick
      self$blast_options  <- blast_options
      self$read_file()
    }, # initialize

    #' @description Read a FASTA file
    #' @param ... further arguments for \code{read_fastafile}
    read_file = function(...){
      self <- read_fastafile(self, ...)
    }, # read_file
    
    #' @description tabulate tetramers, compute PCA and select outliers 
    #' @param verbose logical, if TRUE output some chatter
    tabulate_tetramers = function(verbose = FALSE){
      rawnames <- names(self$rawdata)
      # tabulate
      if (verbose[1]) cat("tabulating windows\n")
      Z <- tabulate_windows(self$rawdata,   
          tWindow = self$parameters[['WINDOW']],
          tStep = self$parameters[['STEP']],
          tWidth = self$parameters[['WIDTH']])
      # differentiate between failed and successful tabulations
      if (verbose[1]) cat("idebtifying failed tabulations\n")
      ok <- sapply(Z, "[[", "status") == 0
      X <- lapply(Z[ok], "[[", "x")
      self$fail <- Z[!ok]
      # condense reverse compliments
      if (verbose) cat("condensing reverse complements\n")
      X <- lapply(X, reduce_tetra)
      # normalize the counts
      if (verbose[1]) cat("normalizing counts\n")
      X <- lapply(X, normalize_tetra)
      if (verbose[1]) cat("binding tabulations into one matrix\n")
      # bind into one matrix
      self$data <- bind_tetra(X)
      # compute PCs
      if (verbose[1]) cat("computing principal components\n")
      self$PC <- pca_tetra(self$data)
      # select outliers
      if (verbose[1]) cat("selecting outliers\n")
      self$outliers <- select_outliers(self)
      # save results
      if (verbose[1]) cat("writing results\n")
      write_tetra(self, verbose = verbose)
    },
    
    #' @description Retrieve the normalized tabulations as a matrix or table
    #' @param form character, either 'matrix' or 'table' (default) 
    #' @return matrix or table of normalized tetramer counts
    get_tabulation = function(form = c("matrix", "table")[2]){
      
      if (tolower(form[1]) == 'table'){
        x <- dplyr::tibble(wname = rownames(X$data)) %>%
          dplyr::bind_cols(dplyr::as_tibble(X$data))
      } else {
        x <- X$data
      }
      x
    },
    
    #' @description Retrieve the principal components as a matrix or table
    # @param form character, either 'matrix' or 'table' (default) 
    #' @param npc numeric, the first npc components are returned
    #' @return matrix or table of principle components
    get_principal_components = function(
      form = c("matrix", "table")[2],
      npc = 8){
      
      if (tolower(form[1]) == 'table'){
        wname <- rownames(X$PC$x)
        cname <- unname(decomposeContigNames(wname)[,1])
        x <- dplyr::tibble(cname, wname) %>%
          dplyr::bind_cols(dplyr::as_tibble(X$PC$x[,seq_len(npc)]))
      } else {
        x <- X$PC$x[,seq_len(npc)]
      }
      x
    },
    
    #' @description Retrieve the principal component rotational loadings as a matrix or table
    # @param form character, either 'matrix' or 'table' (default) 
    # @param npc numeric, the first npc components are returned
    #' @return matrix or table of principle component roational loadings
    get_loadings = function(
      form = c("matrix", "table")[2],
      npc = 8){
      
      if (tolower(form[1]) == 'table'){
        x <- dplyr::tibble(wname = rownames(X$PC$rotation)) %>%
          dplyr::bind_cols(dplyr::as_tibble(X$PC$rotation[,seq_len(npc)]))
      } else {
        x <- X$PC$rotation[,seq_len(npc)]
      }
      x
    },
    
    #' @description Retrieve outlier sequences
    #' @return a table of [cname, wname, PC, pick, sequence] which is an
    #'   augmented version of the outlers table
    get_outlier_sequences = function(){
      self$outliers %>%
        dplyr::mutate(seq = unname(sapply(.data$wname, function(wn) getBlastWindow(self, wn))))
    },
    
    #' @description load blast results
    #' @param filename character, xml filename
    #' @return logical, TRUE if successful
    load_blast = function(filename = file.path(self$out_dir, paste0(self$name,"-outliers.xml"))){
      y <- read_blast_output(filename)
      OK <- TRUE
      if (!is.null(y$info)) {
        self$blast_info <- y$info
      } else {
        OK <- FALSE
      }
      if (!is.null(y$data)) {
        self$blast_data <- y$data
      } else {
        OK <- FALSE
      }
      invisible(OK)
    },
    
    #' @description Run blast on the outlier fasta
    # @param verbose logical, if TRUE output messages for debugging purposes
    run_blast = function(verbose = FALSE){
      if (verbose[1]) cat("calling run_blast_app()\n")
      run_blast_app(self, verbose = verbose)
    }, 
    
    #' @description Run blast on the outlier fasta or simply load existsing blast ouptut
    # @param verbose logical, if TRUE output messages for debugging purposes
    run_or_load_blast = function(verbose = FALSE){
      blastfile <- file.path(self$out_dir, paste0(self$name, "-outliers.xml"))
      if (file.exists(blastfile)){
        if (verbose[1]) cat("loading existing blast data\n")
        self$load_blast(filename = blastfile)
      } else {
        self$run_blast(verbose = TRUE)
      }
    }, 
    
    
    #' @description Run blast on the outlier fasta
    #' @param verbose logical, if TRUE output messages for debugging purposes
    plot = function(addBlast = !is.null(self$blast_data)){
      plot_tetramers(self)
    }  
    
    ) # end of public
  )