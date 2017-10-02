
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
  negIdx <- which( GenomicRanges::start(gr) <= 0 )

  outIdx <- union(outOfBoundIdx, negIdx)

  # save length of left and right out-of-bound lengths
  outGR <- gr[outIdx]

  chromEnds <- GenomeInfoDb::seqlengths(outGR)[
    as.vector(GenomeInfoDb::seqnames(outGR))]

  outDF <- data.frame(
    idx = outIdx,
    left = ifelse(
      GenomicRanges::start(outGR) <= 0,
      abs(GenomicRanges::start(outGR)) + 1,
      0),
    right = ifelse(
      !is.na(chromEnds) &
      GenomicRanges::end(outGR) >
        chromEnds,
      abs(GenomicRanges::end(outGR) - chromEnds),
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


#' Add coverage to regions in \code{\link[GenomicRanges]{GRanges}} object.
#'
#' This function adds a vector of coverage (or any other signal in the input
#' bigWig file) to each range in a \code{\link[GenomicRanges]{GRanges}} object.
#' The coverage is reported for a fixed-sized window around the region center.
#' For regions with negative strand, the coverage vector is reversed.
#'
#' @param gr \code{\link[GenomicRanges]{GRanges}} object with genomic regions
#' @param bwFile File path or connection to BigWig file with coverage to parse
#'   from.
#' @param window the window size around the center of ranges in \code{gr}.
#' @param binSize size of bins to which the coverage values are combined.
#' @param colname name of the new column that is created in \code{gr}.
#'
#' @return \code{\link[GenomicRanges]{GRanges}} as input but with an additional
#'   meta column containing the coverage values for each region.
#' @export
addCovToGR <- function(gr, bwFile, window = 1000, binSize = 1, colname = "cov"){

  # get windows around gr without warnding if ranges extend chromosome borders
  suppressWarnings(
    ancWin <- GenomicRanges::resize(gr, width = window, fix = "center")
  )

  # Because quering coverage result in error for ranges outisde chromosome we
  # need to get intervals defined outside of chromosomes
  outDF <- getOutOfBound(ancWin)

  # trim ranges to fint within chromosomes
  ancWin <- GenomicRanges::trim(ancWin)

  # trim start in case there is no seqinof object
  GenomicRanges::start(ancWin) <- ifelse(
    GenomicRanges::start(ancWin) > 0,
    GenomicRanges::start(ancWin),
    1)

  # get numeric with coverage of each region

  # define query region and trim seqinfo to avoid
  # warning in rtracklayer::import.bw
  selectWin <- ancWin
  GenomeInfoDb::seqlevels(selectWin) <- as.character(
    unique(GenomeInfoDb::seqnames(selectWin))
    )
  GenomeInfoDb::seqinfo(selectWin) <- GenomeInfoDb::Seqinfo(
    GenomeInfoDb::seqlevels(selectWin)
    )
  selection <- rtracklayer::BigWigSelection(selectWin)

  covGR <- rtracklayer::import.bw(
    bwFile,
    selection = selection,
    as = "GRanges",
    seqinfo = seqinfo(ancWin))

  # update covGR with seqinfo to allow subsetting with ancWin
  if ( !any(is.na(GenomeInfoDb::seqlengths(ancWin))) ) {
    GenomeInfoDb::seqlevels(covGR) <- GenomeInfoDb::seqlevels(ancWin)
    GenomeInfoDb::seqinfo(covGR) <- GenomeInfoDb::seqinfo(ancWin)
  }

  covRle <- GenomicRanges::coverage(covGR, weight = covGR$score)

  # get coverage as for all regions
  covAncRle <- covRle[ancWin]
  covAnc <- IRanges::NumericList(covAncRle)

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
    covAnc <- IRanges::NumericList(lapply(covAnc, slideMean, k = binSize))
  }

  # reverse coverage vector for regons on minus strand
  negStrand <- which(as.logical(GenomicRanges::strand(gr) == "-"))
  covAnc[negStrand] <- lapply(covAnc[negStrand], rev)

  # add as additional column to GRanges object
  S4Vectors::mcols(gr)[, colname] <- covAnc

  return(gr)
}
