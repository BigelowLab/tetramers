## tetramers

`tetramers` is an R package detects outlier regions in genomes in terms of tetramer composition, which are indicative of either contamination or gene transfer.  Use is either interactove or via a script included with the package.


### Requirements:

 In addition to [R](https://www.r-project.org/) version 3.0+ and [blast](https://blast.ncbi.nlm.nih.gov) the following packages are required.
  
  + [rlang](https://cran.r-project.org/package=rlang)
  
  + [R6](https://cran.r-project.org/package=R6)
  
  + [dplyr](https://cran.r-project.org/package=dplyr)
   
  + [readr](https://cran.r-project.org/package=readr)
  
  + [yaml](https://cran.r-project.org/package=yaml)
  
  + [blastxml](https://github.com/BigelowLab/blastxml) - installed via github, see below
 
  + [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
 
  + [RColorBrewer](https://cran.r-project.org/package=RColorBrewer)


### Database Requirement

The database must be from 2020 or prior.  
 
### Installation

Choose the library path - see script notes below.  For example...

```r
lib_path = c(
    "/mnt/scgc/scgc_nfs/opt/common/anaconda/2.1.0/lib/R/library",
    "/mnt/scgc/scgc_nfs/opt/common/anaconda3/4.0.0/lib/R/library"
    )
```


Install CRAN packages using `install.packages` . Only do this is you don't have the packages already.

```r
install.packages(c("rlang", "RColorBrewer", "dplyr", "readr", "yaml"), lib = lib_path)
```

Install Bioconductor packages using `bioclite`.  Only do this is you don't have the package already.

```r
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings", lib = lib_path)
``` 

Install github packages using [remotes](https://cran.r-project.org/package=remotes) package.  You must be credentialed to do this step. 

```r
remotes::install_github("BigelowLab/blastxml", lib = lib_path)
remotes::install_github("BigelowLab/tetramers", lib = lib_path)
``` 


### Interactive usage

The package exposes the essentialy steps for teteramer analysis, outlier selection and plotting.

```r
library(tetramers)
path = "/mnt/scgc_nfs/lab/ben/tetramer_test/AH-906-G21_contigs/tetramers"
filename = "AH-906_G21_contigs.fasta"
fasta_file <- file.path(path, filename)
X = Tetramers$new(fasta_file)
X$tabulate_tetramers(verbose = TRUE)
X$run_or_load_blast(verbose = TRUE)
X$plot()
```


# TO DO - update script

### Script location

The script is located in the `scripts` subdirectory for tetramerpca library directory. Library paths within R can be found using `.libPaths()` but here are two examples.

 + `/mnt/scgc/scgc_nfs/opt/common/anaconda/2.1.0/lib/R/library/tetramerpca/scripts/tetramer_pipeline.Rscript`

 + `/mnt/scgc/scgc_nfs/opt/common/anaconda3/4.0.0/lib/R/library/tetramerpca/scripts/tetramer_pipeline.Rscript`
 
 + `/mnt/scgc_nfs/opt/common/R/3.2.2/lib64/R/library/tetramerpca/scripts/tetramer_pipeline.Rscript`
 
Note that the script can also be copied from [github](https://github.com/BigelowLab/tetramers/blob/master/inst/scripts/tetramer_pipeline.Rscript) and placed anywhere convenient.

### Script usage

```
Usage:
tetramer_pipeline.Rscript [--input character] [--output_dir character]
     [--window numeric] [--step numeric] [--blast_cmd character] [--db
     character] [--num_threads character] [--num_alignments character]
     [--evalue character]

Argument details follow
--input (required, type 'character') 
    required input fasta file 
    default:  
--output_dir (type 'character') 
    optional output directory, by default directory same as input 
    default:  
--window (type 'numeric') 
    window width in characters, default is 1600 
    default: 1600 
--step (type 'numeric') 
    window step size in characters, default is 200 
    default: 200 
--blast_cmd (type 'character') 
    name of the bast application possibly with path 
    default: blastn 
--db (type 'character') 
    blast database name possibly with path 
    default: nt 
--num_threads (type 'character') 
    blast option passed through 
    default: 12 
--num_alignments (type 'character') 
    bast option passed through 
    default: 10 
--evalue (type 'character') 
    bast option passed through 
    default: 10 
```

Another example...

```
Rscript --vanilla /mnt/scgc_nfs/opt/common/R/3.2.2/lib64/R/library/tetramerpca/scripts/tetramer_pipeline.Rscript --input something.fasta --output_dir somewhere --num_threads 12 --window 1600 --step 200
```


### What the process does:

 + Generates a matrix of tetramer frequencies (columns) of each contig window (rows), using the window and step sizes. The number of tetramers is 136, reduced from the maximum of 256 to eliminate reverse-complementary tetramers.
  
 + Normalizes the tetramer frequency table as an input to estimate the first eight principal components.

 + For each of the principal components, identifies those contigs that contain "outliers", i.e. windows with either extreme positive or negative values. 

 + Generate blastn output for the selected outliers.

 + Generates four scatter-plots: PC1 vs PC2, PC3 vs PC4, PC5 vs PC6 and PC7 vs PC8. All dots are grey except for those corresponding to contigs that contain outlier windows. Each outlier-containing contig is plotted using separate color and/or symbol, so that they are clearly visible on the plot.  A color/symbol legend is generated for each plot.

 + Output a variety of files where `name` is from the input file.
 
   - Original FASTA sequence, `name.fasta`
   
  - Principal component values, `name-tetramer-PC.csv.gz`
  
  - Principla component loadings, `name-tetramer-loading.csv.gz`
  
  - Listing of contigs where tetramer analysis failed,  `name-tetramer-fail.csv`
  
  - Listing of  `name-tetramer-counts.csv.gz`
  
  - Blast output, `name-outliers.xml`
  
  - Sequence data for outliers (fed to blast), `name-outliers.fasta`
  
  - Table of info about outliers, `name-outliers.csv`
  
  - Metadata, `name-info.yaml`
  
   - PCA plots, `name-WINDOW-STEP-PC.pdf` (where WINDOW and STEP are integers) 
