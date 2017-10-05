
#' Get out of chromosomal bound ranges.
#'
#' @param gr A \code{GRanges} object.
#' @return A \code{data.frame} with rows for each range in \code{gr} that
#'   extends out of chromosomes. The first column holds the index of the range in
#'   \code{gr}, the second the size of the overlap to the left of the chromosome
#'   and the third the size of the overlap to the right of the chromosome.
getOutOfBound <- function(gr){

  # check if extended regions are out of chromosome space
  outOfBoundIdx <- GenomicRanges:::get_out_of_bound_index(gr)
  negIdx <- which( start(gr) <= 0 )

  outIdx <- union(outOfBoundIdx, negIdx)

  # save length of left and right out-of-bound lengths
  outGR <- gr[outIdx]

  chromEnds <- seqlengths(outGR)[
    as.vector(seqnames(outGR))]

  outDF <- data.frame(
    idx = outIdx,
    left = ifelse(
      start(outGR) <= 0,
      abs(start(outGR)) + 1,
      0),
    right = ifelse(
      !is.na(chromEnds) &
      end(outGR) >
        chromEnds,
      abs(end(outGR) - chromEnds),
      0)
  )

  return(outDF)

}


#' Sliding mean over x of intervals of size k
#'
#' Source:
#' \url{http://stats.stackexchange.com/questions/3051/mean-of-a-sliding-window-in-r}
#'
#' @param x numeric vector
#' @param k interval size
#' @return numeric vector of length \code{length(x) / k}.
slideMean <- function(x, k){

  n <- length(x)
  spots <- seq(from = 1, to = n, by = k)
  result <- vector(length = length(spots))

  result <- purrr::map_dbl(seq_along(spots), function(i){
    mean(x[seq(spots[i], (spots[i] + k - 1))], na.rm = TRUE)
  })

  return(result)

}


#'Add coverage to regions in \code{\link[GenomicRanges]{GRanges}} object.
#'
#'This function adds a \code{\link[IRanges]{NumericList}}  of coverage (or any
#'other signal in the input bigWig file) to each range in a
#'\code{\link[GenomicRanges]{GRanges}} object. The coverage is reported for a
#'fixed-sized window around the region center. For regions with negative strand,
#'the coverage vector is reversed. The coverage signal is added as new metadata
#'colum holding a \code{\link[IRanges]{NumericList}} object.
#'This function does not work on Windows. Only NAs are added on Windows.
#'
#'@param gr \code{\link[GenomicRanges]{GRanges}} object with genomic regions
#'@param bwFile File path or connection to BigWig file with coverage to parse
#'  from.
#'@param window Numeric scalar for window size around the center of ranges in
#'  \code{gr}.
#'@param binSize Integer scalar as size of bins to which the coverage values are
#'  combined.
#'@param colname Character as name of the new column that is created in
#'  \code{gr}.
#'
#'@return \code{\link[GenomicRanges]{GRanges}} as input but with an additional
#'  meta column containing the coverage values for each region as
#'  \code{\link[IRanges]{NumericList}}.
#'
#'@examples
#'
#'# use example bigWig file of ChIP-seq signals on human chromosome 22
#'exampleBigWig <- system.file("extdata",
#'"GM12878_Stat1.chr22_1-18000000.bigWig", package = "chromloop")
#'
#'# use example CTCF moitf location on human chromosome 22
#'motifGR <- chromloop::motif.hg19.CTCF.chr22
#'
#'# add ChIP-seq signals to motif regions
#'motifGR <- addCovToGR(motifGR, exampleBigWig)
#'
#'# add ChIP-seq signals as column named "Stat1"
#'motifGR <- addCovToGR(motifGR, exampleBigWig, colname = "Stat1")
#'
#'# add ChIP-seq signals in windows of 500bp around motif centers
#'motifGR <- addCovToGR(motifGR, exampleBigWig, window = 500)
#'
#'# add ChIP-seq signals in bins of 10 bp
#'motifGR <- addCovToGR(motifGR, exampleBigWig, binSize = 10)
#'
#'@import InteractionSet
#'@export
addCovToGR <- function(gr, bwFile, window = 1000, binSize = 1, colname = "cov"){

  # check input arguments
  stopifnot(file.exists(bwFile), length(bwFile) == 1)
  stopifnot(is.numeric(window), length(window) == 1)
  stopifnot(is.numeric(binSize), length(binSize) == 1)
  stopifnot(is.character(colname), length(colname) == 1)

  # check if OS is windows because BioC rtracklayer does not support it BigWig
  # files on windows
  if(.Platform$OS.type == 'windows') {
    warning('Reading of bigWig files is not supported on Winodws. The function addCovToGR() will only add NA.')

    # add NA and return
    mcols(gr)[, colname] <- NA
    return(gr)
  }

  # get windows around gr without warnding if ranges extend chromosome borders
  suppressWarnings(
    ancWin <- resize(gr, width = window, fix = "center")
  )

  # Because quering coverage result in error for ranges outisde chromosome we
  # need to get intervals defined outside of chromosomes
  outDF <- getOutOfBound(ancWin)

  # trim ranges to fint within chromosomes
  ancWin <- trim(ancWin)

  # trim start in case there is no seqinof object
  start(ancWin) <- ifelse(start(ancWin) > 0, start(ancWin), 1)

  # get numeric with coverage of each region

  # define query region and trim seqinfo to avoid
  # warning in rtracklayer::import.bw
  selectWin <- keepSeqlevels(ancWin, unique(seqnames(ancWin)))
  selection <- rtracklayer::BigWigSelection(selectWin)

  covGR <- rtracklayer::import.bw(
    bwFile,
    selection = selection,
    as = "GRanges",
    seqinfo = seqinfo(ancWin))

  # update covGR with seqinfo to allow subsetting with ancWin
  if ( !any(is.na(seqlengths(ancWin))) ) {
    seqlevels(covGR) <- seqlevels(ancWin)
    seqinfo(covGR) <- seqinfo(ancWin)
  }

  covRle <- coverage(covGR, weight = covGR$score)

  # get coverage as for all regions
  covAncRle <- covRle[ancWin]
  covAnc <- NumericList(covAncRle)

  # add NAs for out of bound regions
  covAnc[outDF$idx] <- lapply(seq_along(outDF$idx), function(i){
    covAnc[outDF$idx][[i]] <- c(
      rep(NA, outDF[i, "left"]),
      covAnc[outDF$idx][[i]],
      rep(NA, outDF[i, "right"])
    )
  })

  # combine coverage for bins
  if (binSize > 1) {
    covAnc <- NumericList(lapply(covAnc, slideMean, k = binSize))
  }

  # reverse coverage vector for regons on minus strand
  negStrand <- which(as.logical(strand(gr) == "-"))
  covAnc[negStrand] <- lapply(covAnc[negStrand], rev)

  # add as additional column to GRanges object
  mcols(gr)[, colname] <- covAnc

  return(gr)
}
