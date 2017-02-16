
#  Parse a bigWig file as RLE coverage object
#'
#' The bigWig files contains read counts (or other dense, continouse data)
#' along the genome.
#' The bigWig format is described here:
#' https://genome.ucsc.edu/goldenpath/help/bigWig.html
#'
#' @param inFile Input file path or connection. See \code{con} paramter in
#' \code{\link[rtracklayer]{import}} function.
#' @param seqInfo A \link{\code{seqinfo}} object defining the reference genome.
#' @return An \code{\link{RleList}} object with density values for each
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
#'
addCovToGR <- function(gr, cov, window=1000, bin_size=10){

  # get windows around gr
  ancWin <- GenomicRanges::resize(gr, width=window, fix="center")

  covGRL <- GenomicRanges::ranges(cov)

  ancWinR <- GenomicRanges::restrict(ancWin, covGRL)

  # get numeric with coverage of each region
  covAnc <- IRanges::NumericList(cov[ancWin])


}
