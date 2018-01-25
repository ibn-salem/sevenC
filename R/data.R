#' CTCF motif locations in human genome hg19.
#'
#' A dataset containing the motif hits of the CTCF recognition motif from JASPAR
#' database (MA0139.1, \url{http://jaspar.genereg.net/matrix/MA0139.1/}) in human
#' genome assembly hg19.
#'
#' The dataset was downloaded from JASPAR 2018 motif tracks from the following
#' URL:
# nolint start
#' \url{http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2018/hg19/tsv/MA0139.1.tsv.gz}
# nolint end
#'
#' Motif locations were filtered to contain only motif hists with p-value
#' \eqn{\le 2.5 * 10^-6}.
#' The p-value is the motif hit significance as repoted from the
#' motif scanning alogrithim used during construction of the JASPAR motif
#' tracks. More information on the JASPAR motif track pipeline can be found
#' here:
#' \url{https://github.com/wassermanlab/JASPAR-UCSC-tracks}.
#'
#' @format \code{\link[GenomicRanges]{GRanges}} object with 38774 ranges on
#' positive and negative strand with 1 meta column:
#'
#'   \describe{ \item{score}{
#'   The significance socre of the motif hit, defined as -log_10(p-value).
#'   }
#'   }
#'
#' @source
# nolint start
#' \url{http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2018/hg19/tsv/MA0139.1.tsv.gz}
# nolint end
"motif.hg19.CTCF"

#' CTCF motif locations on chromosome 22 in human genome hg19.
#'
#' A dataset containing the motif hits of the CTCF recognition motif from JASPAR
#' database (MA0139.1) in human genome assembly hg19. Only motifs with a p-value
#' \eqn{\le 2.5 * 10^-6} on chromosome 22 are reported.
#'
#' See '?motif.hg19.CTCF' for a more details and the full data set.
#'
#' @seealso \code{\link{motif.hg19.CTCF}}
"motif.hg19.CTCF.chr22"

#' CTCF motifs on human chromosome 22 with example coverage.
#'
#' This dataset is the same as \code{\link{motif.hg19.CTCF.chr22}} but with one
#' additional metadata colum, named "chip", holding ChIP-seq signals for all
#' motifs in a windows of 1000 bp around the motif center as
#' \code{\link[IRanges]{NumericList}}. The data is from a ChIP-seq experiment
#' for STAT1 in human GM12878 cells. The full bigWig file can be downloaded from
#' ENCODE (Dunham et al. 2012)
# nolint start
#' \url{http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878Stat1StdSig.bigWig}.
# nolint end
#' See \code{\link{motif.hg19.CTCF}} and \code{\link{motif.hg19.CTCF.chr22}} for
#' a more details and the motif data set.
#'
#' @seealso \code{\link{motif.hg19.CTCF}}, \code{\link{motif.hg19.CTCF.chr22}}
#'
"motif.hg19.CTCF.chr22.cov"


#' Default parameters for logistic regression model in sevenC.
#'
#' This dataset contains term names and estimates for logistic regression model
#' to predict chromatin looping interactions. The estimate represent an average
#' of the 10 best performing models out of 124 transcription factor ChIP-seq
#' data sets from ENCODE.
#'
#' Each of 124 transcription factor (TF) ChIP-seq data sets from ENCODE in
#' GM12878 cells were used to train a logistic regression model. All CTCF motifs
#' in \code{\link{motif.hg19.CTCF}} within a distance of 1 Mb were used as
#' candidates. A given pair was labled as true loop interactions, if it has
#' interaction support based on loops from Hi-C in human GM12878 cells from Rao
#' et al. 2014 or ChIA-PET loops from Tang et al. 2015 in the same cell type.
#' The 10 best performing models were selected based on the average area under
#' the precision-recall-curve in 10-fold cross-validation.
#' The parameters were than averaged across the 10 best performig models.
#'
#'
#' @references
#'
#' Suhas S.P. Rao, Miriam H. Huntley, Neva C. Durand, Elena K. Stamenova, Ivan
#' D. Bochkov, James T. Robinson, Adrian L. Sanborn, Ido Machol, Arina D. Omer,
#' Eric S. Lander, Erez Lieberman Aiden, A 3D Map of the Human Genome at
#' Kilobase Resolution Reveals #' Principles of Chromatin Looping, Cell, Volume
#' 159, Issue 7, 18 December 2014, Pages 1665-1680, ISSN 0092-8674,
#' https://doi.org/10.1016/j.cell.2014.11.021.
#'
#' Zhonghui Tang, Oscar Junhong Luo, Xingwang Li, Meizhen Zheng, Jacqueline
#' Jufen Zhu, Przemyslaw Szalaj, Pawel Trzaskoma, Adriana Magalska, Jakub
#' Wlodarczyk, Blazej Ruszczycki, Paul Michalski, Emaly Piecuch, Ping Wang,
#' Danjuan Wang, Simon Zhongyuan Tian, May Penrad-Mobayed, Laurent M. Sachs,
#' Xiaoan Ruan, Chia-Lin Wei, Edison T. Liu, Grzegorz M. Wilczynski, Dariusz
#' Plewczynski, Guoliang Li, Yijun Ruan, CTCF-Mediated Human 3D Genome
#' Architecture Reveals Chromatin Topology for Transcription, Cell, Volume 163,
#' Issue 7, 17 December 2015, Pages 1611-1627, ISSN 0092-8674,
#' https://doi.org/10.1016/j.cell.2015.11.024.
#'
#' @format An object of class \code{data.frame} with 7 rows and 2 columns
#'   holding the term name and estimate.
#'
#'   \describe{
#'
#'   \item{(Intercept)}{The intercept of the logistic regression model.}
#'
#'   \item{dist}{The genomic distance between the centers of motifs in base
#'   pairs (bp).}
#'
#'   \item{strandOrientationdivergent}{Orientation of motif pairs. 1 if
#'   divergent 0 if not.}
#'
#'   \item{strandOrientationforward}{Orientation of motif pairs. 1 if forward 0
#'   if not.}
#'
#'   \item{strandOrientationreverse}{Orientation of motif pairs. 1 if reverse 0
#'   if not.}
#'
#'   \item{score_min}{Minimum of motif hit score between both motifs in pair.
#'   The motif score is defined as -log_10 of the p-value of the motif hit as
#'   reported by JASPAR motif tracks. The unit is -log_10(p) where p is the
#'   p-value of the motif hit.}
#'
#'   \item{cor}{Pearson correlation coefficient of ChIP-seq signals across +/-
#'   500 bp around CTCF motif centers.}
#'
#'   }
#'
#' @seealso \code{\link{cutoffBest10}} and \code{\link{TFspecificModels}}
"modelBest10Avg"

#' TF specific parameters for logistic regression in sevenC
#'
#' sevenC was trained on 124 TF ChIP-seq data sets from ENCODE. Specific
#' parameters are provided in this data set.
#'
#' @format A \code{data.frame} with 868 rows and 7 columns.
#'
#'   \describe{
#'
#'   \item{TF}{TF name used in ChIP-seq experiment.}
#'
#'   \item{file_accession}{File accession ID from ENCODE project}
#'
#'   \item{term}{Model term name. See \code{\link{modelBest10Avg}} for more
#'   detials.}
#'
#'   \item{estimate_mean}{Mean parameter estimate in 10-fold cross-validation}
#'
#'   \item{estimate_median}{Median parameter estimate in 10-fold
#'   cross-validation}
#'
#'   \item{estimate_sd}{Standard deviation of parameter estimate in 10-fold
#'   cross-validation}
#'   }
#'
#' @seealso \code{\link{modelBest10Avg}} and \code{\link{cutoffByTF}}
"TFspecificModels"

#' Optimal cutoff values for logistic regression models.
#'
#' This dataset contains optimal cutoff scores for the response value of logistic
#' regression models. The cutoff is based on optimal F1-scores. A separate model
#' was trained For each of 124 TF ChIP-seq datasets in human GM12878 cells. The
#' model performance was calculated with Hi-C and ChIA-PET interactions using
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
#'@seealso \code{\link{modelBest10Avg}} and \code{\link{TFspecificModels}}
"cutoffByTF"

#' Default optimal cutoff value of logistic regression.
#'
#' This value is the average optimal cutoff value on the 10 best performing TF
#' ChIP-seq data sets. It is used as default cutoff value on the logistic
#' regression response score in \code{\link{predLoops}} function See
#' \code{?'cutoffByTF'} for more details.
#'
#' @seealso \code{\link{cutoffByTF}}, \code{\link{modelBest10Avg}}
"cutoffBest10"
