% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strings.R
\name{decomposeContigNames}
\alias{decomposeContigNames}
\title{Parse the name of a contig}
\usage{
decomposeContigNames(x, sep = "_")
}
\arguments{
\item{x}{character vector of contig names}

\item{sep}{character separator around which the component names are split}
}
\value{
a nrowx2 character array
}
\description{
Given the contig-window names "blahblah_1-5000" return a 2 element
character vector of c("blahblah", "1-5000")
}
