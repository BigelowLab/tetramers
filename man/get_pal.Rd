% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphics.R
\name{get_pal}
\alias{get_pal}
\title{Returns the even-numbered color palette associated with the count provided.}
\usage{
get_pal(x, name = "PuOr")
}
\arguments{
\item{x}{- the number of color steps, if odd then n+1 is returned,
if n > maxN (see \code{RColorBrewer}) the maxN colors are returned with a
warning otherwise n colors are returned}

\item{name}{the name of the ColorBrewer palette to return}
}
\value{
the ColorBrewer palette specified
}
\description{
Returns the even-numbered color palette associated with the count provided.
}
