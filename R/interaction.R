#' Build a \code{\link{GInteractions}} object with all pairs of input
#' \code{\link{GRanges}} within a given distance.
#'
#' Distance is calculated from the center of input regions.
#'
#' @param inGR \code{\link{GRanges}} object of genomic regions. The ranges shuld
#'   be sorted according to chr, strand, and start position. Use
#'   \code{\link[GenomicRanges]{sort}} to sort it.
#' @param maxDist maximal distance in base-pairs between pairs of ranges as
#'   single  numeric value.
#' @return A \code{\link[InteractionSet]{GInteractions}} object with all pairs
#'   within the given distance.
#' @export
getCisPairs <- function(inGR, maxDist=10^6){

  # check that input GRanges object is sorted
  if( !all(inGR == sort(inGR))) stop("Input ranges inGR need to be sorted. Use
                                     sort(inGR) to sort it")

  # get center postions of each input range
  posGR <- GenomicRanges::resize(inGR, width=1, fix="center")

  # calculate overlap all possible gene pairs within maxDist bp
  hits <- GenomicRanges::findOverlaps(posGR,
                      maxgap=maxDist,
                      drop.redundant=TRUE,
                      drop.self=TRUE,
                      ignore.strand=TRUE)

  # build IntractionSet object
  gi <- InteractionSet::GInteractions(
    S4Vectors::queryHits(hits),
    S4Vectors::subjectHits(hits),
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


#' returns indecies of columns with non-zero variance
#'
#' @param dat data.frame or matirx
#' @return column indecies of columns with non-zero variance
noZeroVar <- function(dat) {
  out <- apply(dat, 2, function(x) length(unique(x)))
  which(out > 1)
}

#' Apply a function to pairs of close genomic regions.
#'
#' This function adds a vector with the resulting functin call for each input
#' interaction. Only works for input interaction within the given \code{maxDist}
#' distance
#'
#'
#' @param gi A sorted \code{\link[InteractionSet]{GInteractions}} object.
#' @param datcol a string matching an annotation column in \code{regions(gi)}.
#'   This collumn is assumed to hold the same number of values for each
#'   interaction \code{NumericList}.
#' @param fun A function that takes two numeric vectors as imput to compute a
#'   summary statsitic. Default is \code{\link[stats]{cor}}.
#' @param colname A string that is used as columnname for the new column in
#'   \code{gi}.
#' @param maxDist maximal distance of pairs in bp as numeric. If maxDist=NULL,
#'   the maximal distance is computed from input interactions gi by
#'   \code{max(pairdist(gi))}.
#' @return A \code{\link[InteractionSet]{GInteractions}} similar to \code{gi}
#'   just wiht an additinoal column added.
#' @export
#' @import data.table
applyToCloseGI <- function(gi, datcol, fun=cor, colname="cor", maxDist=NULL){

  # Algorithm
  # (0) define maxDist
  # (1) Group genome in overlapping bins of size 2*maxDist
  # (2) Run pairwise correlation for all ranges in each bin
  # (3) Combine correlations to data.frame with proper id1 and id2 in first columns
  # (4) Query data frame with input pairs
  #   /uses inner_join() from dplyr like in
  #  http://stackoverflow.com/questions/26596305/match-two-data-frames-based-on-multiple-columns

  # check input
  if ( any(is.na(GenomeInfoDb::seqlengths(gi))) ) stop("gi object need seqlengths.")
  if ( !is(gi, "StrictGInteractions") ) stop("gi should be of class StrictGInteractions")

  #-----------------------------------------------------------------------------
  # (0) define maxDist
  #-----------------------------------------------------------------------------
  if( is.null(maxDist) ){
    maxDist <- max(InteractionSet::pairdist(gi))
  }else{
    if(maxDist < max(InteractionSet::pairdist(gi))) stop("maxDist is smaller than maximal distance between interactions in input gi.")
  }

  #-----------------------------------------------------------------------------
  # (1) group ranges in by bins
  #-----------------------------------------------------------------------------

  message("INFO: Prepare Genomic bins...")

  # if (all(!is.na(seqlengths(ancGR))))
  # create GRanges object for entire genome
  genomeGR <- GenomicRanges::GRanges(GenomeInfoDb::seqinfo(gi))

  # tile genoe in overlapping bins of with 2*maxDist
  binGR <- unlist(GenomicRanges::slidingWindows(genomeGR, 2*maxDist, maxDist))

  hits <- GenomicRanges::findOverlaps(binGR, InteractionSet::regions(gi))

  #-----------------------------------------------------------------------------
  # (2) compute pairwise correlatin for all ranges in each bin
  #-----------------------------------------------------------------------------
  message("INFO: compute correlations for each group...")

  covList <- S4Vectors::mcols(InteractionSet::regions(gi))[,datcol]
  datamat <- as.matrix(covList)

  corMatList <- lapply(1:length(binGR), function(i){

    # #DEBUG:
    # message("DEBUG: index i, ", i)

    # get regions in this bin
    regIdx <- S4Vectors::subjectHits(hits)[S4Vectors::queryHits(hits) == i]

    if (length(regIdx) == 1) {
      dat <- cbind(datamat[regIdx,])
    }else{
      dat <- t(datamat[regIdx,])
    }

    # get indices with non-zero variance (they casue warning and NA in cor())
    subIdx <- noZeroVar(dat)

    n = length(subIdx)

    # compute pairwise correlations for all regions in this bin
    if (n != 1) {

      m <- cor(dat[,subIdx])

    }else{

      m <- 1

    }

    # constract data.table object for all pairs
    corDT <- data.table::data.table(
      rep(regIdx[subIdx], n),
      rep(regIdx[subIdx], each = n),
      array(m)
    )

  })

  #-----------------------------------------------------------------------------
  # (3) combine all data.frames
  #-----------------------------------------------------------------------------
  message("INFO: Combine data.tables of pairwise correlations...")
  # corDF <- data.frame(do.call("rbind", corMatList))
  corDT <- data.table::rbindlist(corMatList)

  #-----------------------------------------------------------------------------
  # (4) Query with input pairs
  #-----------------------------------------------------------------------------

  # names(corDF) <- c("id1", "id2", "val")
  names(corDT) <- c("id1", "id2", "val")
  data.table::setkeyv(corDT, cols = c("id1", "id2"))

  # convert gp into data.table and set keys to id1 and id2 columns
  #names(gp)[1:2] <- c("id1", "id2")
  gpDT <- data.table::data.table(
    id1 = InteractionSet::anchors(gi, type = "first", id = TRUE),
    id2 = InteractionSet::anchors(gi, type = "second", id = TRUE),
    key = c("id1", "id2")
  )

  message("INFO: Query correlation for input pairs...")
  matches <- corDT[gpDT, on=c("id1", "id2"), mult="first"]

  #return(matches$val)
  S4Vectors::mcols(gi)[,colname] <- matches$val

  return(gi)
}


#' Add column to \code{\link{GInteractions}} with overlap support.
#'
#' See overlap methods in \code{\link{InteractionSet}} package for more details
#' on the oberlap calculations: \code{?InteractionSet::overlapsAny}
#'
#' @param gi \code{\link{GInteractions}} object
#' @param subject another \code{\link{GInteractions}} object
#' @param colname name of the new annotation columm in \code{gi}.
#' @param ... addtional arguments passed to \code{\link[IRanges]{overlapsAny}}.
#' @return \code{\link{InteractionSet}} \code{gi} as input but with additonal
#'   annotation column \code{colname} indicationg whether each interaction
#'   is supported by \code{subject} or not.
#' @export
addInteractionSupport <- function(gi, subject, colname = "loop", ...){

  ol <- IRanges::overlapsAny(gi, subject, ...)

  S4Vectors::mcols(gi)[,colname] <- factor(ol,
                                          c(FALSE, TRUE),
                                          c("No loop", "Loop"))

  return(gi)
}


#' Add combination of anchor strand orientation.
#'
#' Each anchor region has a strand that is \code{'+'} or \code{'-'}. Threfore
#' the each interaction between two regions has one of the following strand
#' combinations: "forward", "reverse", "convergent", or "divergent". Unstranded
#' ragnes, indicated by (\code{*}), are treated as positive strand.
#' @param gi \code{\link{GInteractions}}
#' @param colname name of the new colum that is created in \code{gi}.
#'
#' @return The same \code{\link{GInteractions}} as \code{gi} but with an
#'   additonal column indicating the four possible combinations of strands
#'   "forward", "reverse", "convergent", or "divergent".
#' @export
addStrandCombination <- function(gi, colname = "strandOrientation"){

  anc1 <- InteractionSet::anchors(gi, type = "first", id = TRUE)
  anc2 <- InteractionSet::anchors(gi, type = "second", id = TRUE)

  # because interactions in GInteraction are sorted that + strand comes first
  # we need to make sure that the more upstream anchor region is taken as first
  ancStart <- GenomicRanges::start(InteractionSet::regions(gi))
  anc1First <- ancStart[anc1] <= ancStart[anc2]

  first <- ifelse(anc1First, anc1, anc2)
  second <- ifelse(anc1First, anc2, anc1)

  ancStrand <- GenomicRanges::strand(InteractionSet::regions(gi))

  sameStrand <- ancStrand[first] == ancStrand[second]
  firstPos <- as.character(ancStrand[first]) %in% c("+", "*")

  S4Vectors::mcols(gi)[sameStrand & firstPos, colname] <- "forward"
  S4Vectors::mcols(gi)[sameStrand & !firstPos, colname] <- "reverse"
  S4Vectors::mcols(gi)[!sameStrand & firstPos, colname] <- "convergent"
  S4Vectors::mcols(gi)[!sameStrand & !firstPos, colname] <- "divergent"

  return(gi)
}


#' Add motif score of anchors.
#'
#' If each anchor region (motif) has a score as annotation column, this function
#' adds two new columns named "score_1" and "score_2" with the scores of the
#' first and the second anchor region, respectively. Additonally a column named
#' "score_min" is added with holds for each interaction the mimium of "score_1"
#' and "score_2".
#' @param gi \code{\link{GInteractions}}
#' @param colname name of the new colum that is created in \code{gi}.
#' @return The same \code{\link{GInteractions}} as \code{gi} but with three
#'   additonal annotation columns.
#' @export
addMotifScore <- function(gi, colname = "score"){

  anc1 <- InteractionSet::anchors(gi, type = "first", id = TRUE)
  anc2 <- InteractionSet::anchors(gi, type = "second", id = TRUE)

  ancScore <- S4Vectors::mcols(InteractionSet::regions(gi))[, colname]

  S4Vectors::mcols(gi)[, "score_1"] <- ancScore[anc1]
  S4Vectors::mcols(gi)[, "score_2"] <- ancScore[anc2]
  S4Vectors::mcols(gi)[, "score_min"] <- apply(cbind(ancScore[anc1], ancScore[anc2]), 1, min)

  return(gi)
}


#' Get genomic ranges spanned by (intracrhomosomal) interactions.
#'
#' Note, this method assumes are intra-chromosomal (e.g. only interactions
#' between regions from the same chromosome).
#'
#' @param gi \code{\link{GInteractions}} or \code{\link{InteractionSet}} object
#' @return A \code{\link{GRanges}} object with ranges spanning the
#'   interactions in \code{gi}.
#'
#' @export
interactionRange <- function(gi){

  # assume only intra-chromosomal interactions
  stopifnot(all(InteractionSet::intrachr(gi)))

  # make first anchor <= second anchor
  gi <- as(gi, "GInteractions")
  strand(regions(gi)) <- "*"
  gi <- swapAnchors(gi, mode = "order")

  # get chroms from first anchor
  chr <- seqnames(anchors(gi, "first"))

  # take all coordinates to get smalles as start and largest as end
  coords <- list(
    start(anchors(gi, "first")),
    end(anchors(gi, "first")),
    start(anchors(gi, "second")),
    end(anchors(gi, "second"))
  )
  start <- do.call(pmin, coords)
  end <- do.call(pmax, coords)

  # create GenomicRanges object
  gr <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges(start, end),
    strand = "*",
    seqinfo = seqinfo(gi)
  )

  # add mcols
  mcols(gr) <- mcols(gi)

  return(gr)
}

#' Extend interaction anchors with specific term for inner and outer site
#'
#' @param gi \code{\link{GInteractions}} or \code{\link{InteractionSet}} object
#' @param inner,outer A scalar, non-negative, of the extension sizes for anchor
#'   ends outside and inside interaction loops.
#' @return gi \code{\link{GInteractions}} object with extened regions.
extendAnchors <- function(gi, inner, outer){

  # turn gi into GIntreactions
  gi <- as(gi, "GInteractions")
  strand(regions(gi)) <- "*"
  gi <- swapAnchors(gi, mode = "order")

  # extend ranges of anchors
  anc1 <- anchors(gi, "first")
  start(anc1) = start(anc1) - outer
  # extend end coordinate of fist anchor but only until start of second
  end(anc1) = pmin(
    end(anc1) + inner,
    pmax(start(anchors(gi, "second")) - 1, end(anc1))
  )

  anc2 <- anchors(gi, "second")
  # extend start coordinate of second anchor but only until end of first
  start(anc2) = pmax(
    start(anc2) - inner,
    pmin(end(anchors(gi, "first")) + 1, start(anc2))
  )
  end(anc2) = end(anc2) + outer

  extendedGI <- InteractionSet::GInteractions(anc1, anc2)
  mcols(extendedGI) <- mcols(gi)

  return(extendedGI)
}

