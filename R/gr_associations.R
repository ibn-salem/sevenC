
#' Associate two sets of regions together by interactinos.
#'
#' @param gr1 A GRanges object.
#' @param gr2 A GRanges object.
#' @param gi A \code{\link[InteractionSet]{GInteractions}} or
#'   \code{\link[InteractionSet]{InteractionSet}} object.
#' @param ignore.strand Logical. See \code{\link[InteractionSet]{linkOverlaps}}.
#' @param ... Other arguments passed to \code{\link[GenomicRanges]{findOverlaps}}
#'   and \code{\link[InteractionSet]{linkOverlaps}}
#'
#' @return A dataframe of integer indices indicating which elements of
#'   \code{gr1} link which elements of \code{gr1}.
#'
#' @export
linkRegions <- function(gr1, gr2, gi, ignore.strand = TRUE, ...){

  # get direct overlapping elements
  ol <- GenomicRanges::findOverlaps(gr1, gr2, ignore.strand = ignore.strand, ...)

  # link regions by interactions
  links <- InteractionSet::linkOverlaps(gi, gr1, gr2, ignore.strand = ignore.strand, ...)

  # take only gr1 and gr2 indices and rename to match direct overlap
  olDF <- as.data.frame(ol)
  names(olDF) <- c("gr1", "gr2")
  olDF[, "gi"] <- NA

  linksDF <- as.data.frame(links)
  names(linksDF) <-c("gi", "gr1", "gr2")

  # combine direct overlaps and links via interactions
  hits <- rbind(
    olDF,
    linksDF
  )

  return(hits)

}

