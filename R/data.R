#' CTCF motif locations in human genome hg19.
#'
#' A dataset containing the motif hits of the CTCF recogintion motif from JASPAR
#' database (MA0139.1) in human genome assembly hg19. Only motifs with a p-value
#' \eqn{\le 10^-5} are reported.
#'
#' The dataset was created by using the tool RSTA matrix-scan
#' (\url{http://www.rsat.eu/}) to get all motif hits in the hg19 genome sequence.
#' The PFM motif model for CTCF was downloaded from JASPAR database using this
#' url:
#' \url{http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/individual/MA0139.1.pfm}
#' The following parameters where used.
#'
#' \describe{
#'   \item{-quick}{}
#'   \item{origin start}{}
#'   \item{-bginput}{}
#'   \item{-markov 1}{}
#'   \item{-2str}{}
#'   \item{-uth pval 0.0001}{}
#'   \item{-return pvalue}{}
#' }
#' More information on each parameter can be found on the matrix-scan manual
#' page: \url{http://pedagogix-tagc.univ-mrs.fr/rsat/help.matrix-scan.html}
#'
#' The dataset was than filtered to contain only motif hists with
#' p-value \eqn{\le 10^-5}.
#' Any other sentence here.
#'
#' @format GRanges object with 132739 ranges on
#'   positive and negative strnad indicated and 5 meta columns:
#'
#' \describe{
#'   \item{sequence}{genomic DNA sequnece}
#'   \item{weight}{
#'     The program scans the input sequences with a position-specific scoring
#'     matrix (PSSM) by selecting, at each position, a sequence segment (S) of
#'     the same length as the matrix, and assigning a score to this segment.
#'
#'     The segment score (weight) is calculated according to the theory
#'     developed by Jerry Hertz and Gary Stormo (1999), but with the
#'     capability to use Markov chain-based background models, as proposed
#'     by Thijs et al. (2001).
#'
#'     The weight of a sequence segment (Ws) is calculated as the log-ratio
#'     between two probabilities:
#'     \deqn{ Ws = log(P(S|M)/P(S|B)) }
#'     where
#'     \deqn{proba_M = P(S|M)}
#'     is the probability to generate the sequence segment given the matrix.
#'     \deqn{proba_B = P(S|B)}
#'     is the probability to generate the sequence segment given the background
#'     model.}
#'   \item{Pval}{The site-wise P-value estimates the significance of the weight
#'    associated to each site. It is computed from the matrix, according to the
#'     probabilities described in Staden (1989), Bailey (1998).}
#'  \item{ln_Pval}{The natural logarithm of the p-value.}
#'  \item{ln_Pval}{The significance, defined as sig = -log_10(P-value)}
#' }
#' @source \url{http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/individual/MA0139.1.pfm} and \url{http://www.rsat.eu/}
"motif.hg19.CTCF"

