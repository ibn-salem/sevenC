#' CTCF motif locations in human genome hg19.
#'
#' A dataset containing the motif hits of the CTCF recognition motif from JASPAR
#' database (MA0139.1) in human genome assembly hg19.
#'
#' The dataset was created by using the tool RSAT matrix-scan
#' (\url{http://www.rsat.eu/}) to get all motif hits in the hg19 genome
#' sequence. The PFM motif model for CTCF was downloaded from JASPAR database
#' using this url:
# nolint start
#' \url{http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/individual/MA0139.1.pfm}
# nolint end
#' The following parameters where used.
#'
#' \describe{ \item{-quick}{} \item{origin start}{} \item{-bginput}{}
#' \item{-markov 1}{} \item{-2str}{} \item{-uth pval 0.0001}{} \item{-return
#' pvalue}{} } More information on each parameter can be found on the
#' matrix-scan manual page:
#' \url{http://pedagogix-tagc.univ-mrs.fr/rsat/help.matrix-scan.html}
#'
#' The dataset was than filtered to contain only motif hists with p-value
#' \eqn{\le 10^-6}.
#'
#' @format GRanges object with 132739 ranges on positive and negative strand
#'   indicated and 5 meta columns:
#'
#'   \describe{ \item{sequence}{genomic DNA sequence}
#'
#'   \item{weight}{ The program scans the input sequences with a
#'   position-specific scoring matrix (PSSM) by selecting, at each position, a
#'   sequence segment (S) of the same length as the matrix, and assigning a
#'   score to this segment.
#'
#'   The segment score (weight) is calculated according to the theory developed
#'   by Jerry Hertz and Gary Stormo (1999), but with the capability to use
#'   Markov chain-based background models, as proposed by Thijs et al. (2001).
#'
#'   The weight of a sequence segment (Ws) is calculated as the log-ratio
#'   between two probabilities: \deqn{ Ws = log(P(S|M)/P(S|B)) } where
#'   \deqn{proba_M = P(S|M)} is the probability to generate the sequence segment
#'   given the matrix. \deqn{proba_B = P(S|B)} is the probability to generate
#'   the sequence segment given the background model.}
#'
#'   \item{Pval}{The site-wise P-value estimates the significance of the weight
#'   associated to each site. It is computed from the matrix, according to the
#'   probabilities described in Staden (1989), Bailey (1998).}
#'
#'   \item{ln_Pval}{The natural logarithm of the p-value.}
#'
#'   \item{sig}{The significance, defined as sig = -log_10(P-value)} }
#' @source
# nolint start
#' \url{http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/individual/MA0139.1.pfm}
# nolint end
#' and \url{http://www.rsat.eu/}
"motif.hg19.CTCF"

#' CTCF motif locations on chromosome 22 in human genome hg19.
#'
#' A dataset containing the motif hits of the CTCF recognition motif from JASPAR
#' database (MA0139.1) in human genome assembly hg19. Only motifs with a p-value
#' \eqn{\le 10^-6} on chromosome 22 are reported.
#'
#' See '?motif.hg19.CTCF' for a more details and the full data set.
#'
"motif.hg19.CTCF.chr22"

#' CTCF motifs on human chromosome 22 with example coverage.
#'
#' This dataset is the same as \code{\link{motif.hg19.CTCF.chr22}} but with one
#' additional metadata colum, named "chip", holding ChIP-seq signals for all
#' motifs in a windows of 1000 bp around the motif center as
#' \code{\link[IRanges]{NumericList}}. The data is from a ChIP-seq experiment
#' for STAT1 in human GM12878 cells. The full bigWig file can be downloaded from
#' ENCODE (Dunham et al. 2012)
#' [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878Stat1StdSig.bigWig).
#'
#' See \code{\link{motif.hg19.CTCF}} and \code{\link{motif.hg19.CTCF.chr22}} for
#' a more details and the motif data set.
#'
"motif.hg19.CTCF.chr22.cov"


#' Default parameters for logistic regression model.
#'
#' This dataset contain term names and estimates for logistic regression model
#' to predict chromatin looping interactions. The estimate represent an average
#' of the 10 best performing models out of 121 transcription factor ChIP-seq
#' data sets from ENCODE.
#'
#' Each of 121 transcriptin factor (TF) ChIP-seq data sets from ENCODE in
#' GM12878 cells were used to train a logistic regression model. All CTCF motifs
#' in \code{\link{motif.hg19.CTCF}} within a distance of 1 Mb were used as
#' candidates. A given pair was labled as true loop interactions, if it has
#' interaction support based on Hi-C lops in GM12878 from Rao et al. 2014 or
#' ChIA-PET loops from Tang et al. 2015 in GM12878 cells. The 10 best performing
#' models were selected based on the average area under the
#' precision-recall-curve in 10-fold cross-validation.
#'
#' @format An object of class \code{data.frame} with 7 rows and 2 columns:
#'
#'   \describe{
#'
#'   \item{term}{Term name}
#'
#'   \item{estimate}{Parameter estimate}
#'
#'   }
#'
"modelBest10Avg"


#' Optimal cutoff values for logistic regression models.
#'
#' This dataset contains opimal cutoff scores for the response value of logistic
#' regression models. The cutoff is based on optimal F1-scores. A separate model
#' was trained For each of 121 TF ChIP-seq datasets in human GM12878 cells. The
#' model performance were calculated with Hi-C and ChIA-PET interactions using
#' 10-fold cross-validation.
#'
#' @format An object of class \code{tbl_df} with 121 rows and 3 columns:
#'
#'   \describe{
#'
#'   \item{TF}{Transcription factor name}
#'
#'   \item{max_cutoff}{The optimal cutoff on the logistic regression response
#'   value}
#'
#'   \item{max_f1}{The optimal f1-score associated to the \code{max_cutoff}
#'   value}
#'
#'   }
#'
"cutoffByTF"

#' Default optimal cutoff value of logistic regression.
#'
#' This value is the average optimal cutoff value on the 10 best performing TF
#' ChIP-seq data sets. It is used as default cutoff value on the logistic
#' regression response score in \code{\link{predLoops}} fucntion. See
#' \code{?'cutoffByTF'} for more details.
"cutoffBest10"
