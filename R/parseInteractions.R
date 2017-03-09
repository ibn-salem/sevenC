

#' Parse chromatin loops from Rao et al 2014 as an \code{\link{InteractionSet}}
#' object.
#'
#' @param infile input file with loops
#' @param ... additional arguments, that will be passed to
#'   \code{\link{GRanges()}} functions.
#' @return An code{\link{InteractionSet}} with loops from input file.
#' @export
parseLoopsRao <- function(inFile, ...){

  # parse IMR90 domains from Rao et al 2014:
  raoDF = utils::read.delim(inFile)

  # get chromsome column
  chr = paste0("chr", raoDF$chr1)

  # substract 1 bp from down coordinate to have inclusive interval ranges
  upAnchor = GenomicRanges::GRanges(
      chr,
      IRanges::IRanges(raoDF[,"x1"], raoDF[,"x2"]-1),
      ...)

  downAnchor = GenomicRanges::GRanges(
    chr,
    IRanges::IRanges(raoDF[,"y1"], raoDF[,"y2"]-1),
    ...)

  # build GInteractions
  gi <- InteractionSet::GInteractions(upAnchor, downAnchor, mode="strict")

  # define and add additional annotations to each interaction
  annotationCols <- setdiff(names(raoDF),
                            c("chr1", "x1", "x2", "chr2", "y1", "y2"))

  S4Vectors::mcols(gi) <- raoDF[, annotationCols]

  return(gi)
}
