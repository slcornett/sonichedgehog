% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcm_functions.R
\name{phyloAA}
\alias{phyloAA}
\title{phyloAA}
\usage{
phyloAA(fas_msa)
}
\arguments{
\item{msa}{formatted multiple alignment data}
}
\description{
This function conducts the multiple alignment and distance calculations for phylogenetic data
}
\examples{
using the built-in Shh_fasta data of sonic hedgehog signal molecule orthologs
file <- "link to raw data on git hub"
pcms_AlignAA(file) # the output of this will be the sequences aligned by the MUSCLE algorithm
}
\keyword{distance}
\keyword{msa,}
