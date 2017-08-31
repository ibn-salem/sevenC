
#' Associate two sets of regions together by interactions with ajusted anchors
#'
#' @param gr1 A GRanges object.
#' @param gr2 A GRanges object.
#' @param gi A \code{\link[InteractionSet]{GInteractions}} or
#'   \code{\link[InteractionSet]{InteractionSet}} object.
#' @param outer_maxgap A scalar, non-negative, maxial allowed distance outside
#'   interaction loops to be considered as overlap.
#' @param inner_maxgap A scalar, non-negative, maxial allowed distance inside
#'   interaction loops to be considered as overlap.
#' @param ... Other arguments passed to
#'   \code{\link[GenomicRanges]{findOverlaps}} and
#'   \code{\link[InteractionSet]{linkOverlaps}}
#'
#' @return A dataframe of integer indices indicating which elements of
#'   \code{gr1} link which elements of \code{gr1}.
#'
#' @export
linkRegions <- function(gr1, gr2, gi, outer_maxgap = 0, inner_maxgap = 0, ...){

  # build new GInteraction object with extened anchor regions
  gi <- extendAnchors(gi, outer = outer_maxgap, inner = inner_maxgap)

  # get direct overlapping elements
  ol <- GenomicRanges::findOverlaps(gr1, gr2, ...)

  # link regions by interactions
  links <- InteractionSet::linkOverlaps(gi, gr1, gr2, ...)

  # take only gr1 and gr2 indices and rename to match direct overlap
  olDF <- as.data.frame(ol)
  names(olDF) <- c("gr1", "gr2")
  olDF[, "gi"] <- NA

  linksDF <- data.frame(links)
  names(linksDF) <- c("gi", "gr1", "gr2")

  # combine direct overlaps and links via interactions
  hits <- rbind(
    olDF,
    linksDF
  )

  return(hits)
}


#' Associate two sets of regions together when they are within a loop.
#'
#' @param gr1 A GRanges object.
#' @param gr2 A GRanges object.
#' @param gi A \code{\link[InteractionSet]{GInteractions}} or
#'   \code{\link[InteractionSet]{InteractionSet}} object.
#' @param ... Other arguments passed to \code{\link[GenomicRanges]{findOverlaps}}
#'   and \code{\link[InteractionSet]{linkOverlaps}}
#'
#' @return A dataframe of integer indices indicating which elements of
#'   \code{gr1} link which elements of \code{gr1}.
#'
#' @importFrom dplyr select left_join rename
#' @importFrom magrittr '%>%'
#' @export
linkRegionsInLoops <- function(gr1, gr2, gi, ...){

  # get direct overlapping elements
  ol <- GenomicRanges::findOverlaps(gr1, gr2, ...)

  # get range of all looping interactions
  loopGR <- interactionRange(gi)

  gr1Hits <- data.frame(
    GenomicRanges::findOverlaps(gr1, loopGR)
  ) %>%
    dplyr::rename(gr1 = queryHits)

  gr2Hits <- data.frame(
    GenomicRanges::findOverlaps(gr2, loopGR)
  ) %>%
    dplyr::rename(gr2 = queryHits)

  # link reginos by
  linksDF <- dplyr::left_join(gr1Hits, gr2Hits, by = "subjectHits") %>%
    dplyr::select(gi = subjectHits, gr1, gr2)


  # take only gr1 and gr2 indices and rename to match direct overlap
  olDF <- data.frame(ol)
  names(olDF) <- c("gr1", "gr2")
  olDF[, "gi"] <- NA

  # combine direct overlaps and links via interactions
  hits <- rbind(
    olDF,
    linksDF
  )

  return(hits)
}
