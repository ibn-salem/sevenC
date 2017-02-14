#' Build a \code{\link{GInteractions}} object with all pairs of input
#' \code{\link{GRanges}} within a given distance.
#'
#'
#' @param inGR \code{\link{GRanges}} object of genomic regions. The ranges shuld
#'   be sorted according to chr, strand, and start position. Use
#'   \code{\link[GenomicRanges]{sort()}} to sort it.
#' @param maxDist maximal distance in base-pairs between pairs of ranges as
#'   single  numeric value.
#' @return A \code{\link{GInteractions} object with all pairs within the given
#'   distance.
#' @export
getCisPairs <- function(inGR, maxDist=10^6){

  # check that input GRanges object is sorted
  if( !all(inGR == sort(inGR))) stop("Input ranges inGR need to be sorted. Use
                                     sort(inGR) to sort it")

  # calculate overlap all possible gene pairs within maxDist bp
  hits = GenomicRanges::findOverlaps(
      inGR,
      inGR,
      maxgap=maxDist,
      ignore.strand=TRUE)

  # flag pairs with same range
  sameReg <- S4Vectors::queryHits(hits) == S4Vectors::subjectHits(hits)

  gi <- InteractionSet::GInteractions(
    S4Vectors::queryHits(hits)[!sameReg],
    S4Vectors::subjectHits(hits)[!sameReg],
    inGR,
    mode="strict"
  )

  # sort gi
  gi <- BiocGenerics::sort(gi)

  # add distance
  gi$dist <- InteractionSet::pairdist(gi, type="mid")

  # remove pairs with distance >= maxDist
  # this is essential in case of non-zero length ranges in inGR
  gi <- gi[gi$dist <= maxDist]

  return(gi)
}
