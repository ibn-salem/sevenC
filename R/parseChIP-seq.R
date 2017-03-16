
#'  Parse a bigWig file as RleList coverage object
#'
#' The bigWig files contains read counts (or other dense, continouse data)
#' along the genome.
#' The bigWig format is described here:
#' \url{https://genome.ucsc.edu/goldenpath/help/bigWig.html}
#'
#' @param inFile Input file path or connection. See \code{con} paramter in
#' \code{\link[rtracklayer]{import}} function.
#' @param seqInfo A \code{\link{seqinfo}} object defining the reference genome.
#' @return An \code{\link[IRanges]{RleList}} object with density values for each
#' position in the genome.
#' @export
parseBigWigToRle <- function(inFile, seqInfo, format="bigWig"){

  # parse file as GRange object
  bwGR = rtracklayer::import(inFile, format=format)

  # set seqinfo object
  GenomeInfoDb::seqlevels(bwGR) <- GenomeInfoDb::seqlevels(seqInfo)
  GenomeInfoDb::seqinfo(bwGR) <- seqInfo

  # compute coverage by using the score as weight
  bwCov = GenomicRanges::coverage(bwGR, weight="score")
  return(bwCov)
}


#' Add coverage to regions in \code{\link[GenomicRanges]{GRanges}} object.
#'
#' This function adds a coverage vector to each range in a genomic ranges
#' object. The coverage is reported for a fixed-sized window around the region
#' center and is reversed in case of negative strand of the region.
#'
#' @param gr \code{\link[GenomicRanges]{GRanges}} object with genomic regions
#' @param cov \code{\link{RleList}} object with coverage for each chromosome. Such
#'   an object is returned from \code{\link{parseBigWigToRle}} function.
#' @param window the window size arund the center of ranges in \code{gr}.
#' @param bin_size size of bins to which the coverage values are combined.
#' @return \code{\link[GenomicRanges]{GRanges}} as input but with an additional meta column containing the coverage values for each region.
#' @export
addCovToGR <- function(gr, cov, window=1000, bin_size=10, colname="cov"){

  # get windows around gr
  suppressWarnings(ancWin <- GenomicRanges::resize(gr, width=window, fix="center"))

  # check if extended regions are out of chromsome space
  outOfBoundIdx <- GenomicRanges:::get_out_of_bound_index(ancWin)

  # save length of left and right out-of-bound lengths
  outGR <- ancWin[outOfBoundIdx]

  outDF <- data.frame(
    idx = outOfBoundIdx,
    left = ifelse(GenomicRanges::start(outGR) <= 0,
                  abs(GenomicRanges::start(outGR)) + 1,
                  0),
    right = ifelse(GenomicRanges::end(outGR) > GenomeInfoDb::seqlengths(outGR),
                   abs(GenomicRanges::end(outGR) - GenomeInfoDb::seqlengths(outGR)),
                   0)
  )

  # if (length(outOfBoundIdx) > 0) {
  #   stop("Windows around regions extend out of chromosomal bounds.")
  # }

  # # check taht window reg does not extend range of coverage data
  # covRange <- sapply(cov, length)
  # winRange <- range(range(ancWin))

  # trim ranges to fint within chromosomes
  ancWin <- GenomicRanges::trim(ancWin)

  # get numeric with coverage of each region
  covAnc <- IRanges::NumericList(cov[ancWin])

  # add NAs for out of bound regions
  covAnc[outDF$idx] <- lapply(seq_along(outDF$idx), function(i){
    covAnc[outDF$idx][[i]] <- c(
      rep(NA, outDF[i, "left"]),
      covAnc[outDF$idx][[i]],
      rep(NA, outDF[i, "right"])
    )
  })

  # reverse coverage vector for regons on minus strand
  # TODO

  # add as additional column to GRanges object
  S4Vectors::mcols(gr)[,colname] <- covAnc

  return(gr)

}
