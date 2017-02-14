
#  Parse a bigWig file as RLE coverage object
#'
#' The bigWig files contains read counts (or other dense, continouse data)
#' along the genome.
#' The bigWig format is described here:
#' https://genome.ucsc.edu/goldenpath/help/bigWig.html
#'
#' @param inFile Input file path or connection. See \code{con} paramter in
#' \code{\link[rtracklayer]{import}} function.
#' @param seqInfo A \cod\link{seqinfo} object defining the reference genome.
#' @return An RLE object with density values for each position in the genome.
#' @export
parseBigWigToRle <- function(inFile, seqInfo, format="bigWig"){

  # parse file as GRange object
  bwGR = rtracklayer::import(inFile, format=format)

  # set seqinfo object
  GenomeInfoDB::seqlevels(bwGR) <- GenomeInfoDB::seqlevels(seqInfo)
  GenomeInfoDB::seqinfo(bwGR) <- seqInfo

  # compute coverage by using the score as weight
  bwCov = GenomicRanges::coverage(bwGR, weight="score")
  return(bwCov)
}
