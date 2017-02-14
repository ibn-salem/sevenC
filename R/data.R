#' CTCF motif locations in human genome hg19.
#'
#' A dataset containing the motif hits of the CTCF recogintion motif from JASPAR
#' database (MA0139.1) in human genome assembly hg19. Only motifs with a p-value
#' smaller or equal to 0.0001 are reported.
#'
#' The dataset was created by using the tool RSTA matrix-scan
#' (\url{http://www.rsat.eu/}) to get all motif hits in the hg19 genome sequence.
#' The PFM motif model for CTCF was downloaded from JASPAR database using this
#' url:
#' \url{http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/individual/MA0139.1.pfm}
#' The following parameters where used.
#'
#' \describe{
#'   \item{-quci}{}
#'   \item{origin start}{}
#'   \item{-bginput}{}
#'   \item{-markov 1}{}
#'   \item{-2str}{}
#'   \item{-uth pval 0.0001}{}
#'   \item{-return pvalue}{}
#' }
#'
#' @format \link[GenomicRanges]{\code{GRanges}} object with 132739 ranges on
#'   positive and negative strnad indicated and 5 meta columns:
#' \describe{
#'   \item{sequence}{genomic DNA sequnece}
#'   \item{weight}{weight of the diamond, in carats}
#' }
#' @source \url{http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/individual/MA0139.1.pfm} and \url{http://www.rsat.eu/}
"motif.hg19.CTCF"

