% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pfm_comparison.R
\name{motifSimilarity}
\alias{motifSimilarity}
\title{Calculate similarities between pairs of motifs.}
\usage{
motifSimilarity(
  x,
  y = NULL,
  method = c("R", "HOMER"),
  homerfile = findHomer("compareMotifs.pl"),
  homerOutfile = NULL,
  Ncpu = 1L,
  verbose = TRUE
)
}
\arguments{
\item{x}{Either a \code{\link[TFBSTools]{PFMatrixList}}, or a character
scalar with a file containing motifs in HOMER format (used directly
\code{method = "HOMER"}, loaded into a
\code{\link[TFBSTools]{PFMatrixList}} by \code{\link{homerToPFMatrixList}}
for \code{method = "R"}).}

\item{y}{Either a \code{\link[TFBSTools]{PFMatrixList}} or \code{NULL}
(default). If \code{y = NULL}, then similarities will be calucalted for all
pairs of motifs within \code{x}. Otherwise, \code{method} must be
\code{"R"} and similarities will be calculated between any motif from
\code{x} to any motif from \code{y}.}

\item{method}{A character scalar specifying the method for similarity
calculations. Either \code{"R"} (pure R implementation) or \code{"HOMER"}
(will call the \code{compareMotifs.pl} script from HOMER). Results are
identical (appart from rounding errors), and the R implementation is
usually faster and can be parallelized (\code{Ncpu} argument).}

\item{homerfile}{Path to the HOMER script \code{compareMotifs.pl} (only used
for \code{method = "HOMER"}.}

\item{homerOutfile}{A character scalar giving the file to save the similarity
scores (only for \code{metho = "HOMER"}). If \code{NULL}, scores will be
stored into a temporary file.}

\item{Ncpu}{The number of CPU cores to use when calculating similarities.
This uses \code{\link[parallel]{mclapply}} and only works for \code{method
= "R"}.}

\item{verbose}{A logical scalar. If \code{TRUE}, report on progress.}
}
\value{
A matrix of Pearson's correlation coefficients for each pair of
  motifs.
}
\description{
For each pair of motifs, calculate the similarity defined as the
  maximal Pearson's correlation coefficient between base frequencies over all
  possible shifts (relative positions of the two matrices with at least one
  overlapping position). If necessary matrices are padded on the sides with
  background base frequencies (assuming all bases to have a frequency of
  0.25) to enable comparison of all positions in both matrices.
}
\seealso{
\code{\link[parallel]{mclapply}} for how parallelization is done,
  documentation of HOMER's \code{compareMotifs.pl} for details on
  \code{method = "HOMER"}.
}