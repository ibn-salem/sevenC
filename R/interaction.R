#' Build a \code{\link{GInteractions}} object with all pairs of input
#' \code{\link[GenomicRanges]{GRanges}} within a given distance.
#'
#' Distance is calculated from the center of input regions.
#'
#' @param inGR \code{\link[GenomicRanges]{GRanges}} object of genomic regions. The ranges shuld
#'   be sorted according to chr, strand, and start position. Use
#'   \code{\link[BiocGenerics]{sort}} to sort it.
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
                      maxgap = maxDist,
                      drop.redundant = TRUE,
                      drop.self = TRUE,
                      ignore.strand = TRUE)

  # build IntractionSet object
  gi <- InteractionSet::GInteractions(
    S4Vectors::queryHits(hits),
    S4Vectors::subjectHits(hits),
    inGR,
    mode = "strict"
  )

  # sort gi
  gi <- BiocGenerics::sort(gi)


  # add distance
  gi$dist <- InteractionSet::pairdist(gi, type = "mid")

  # add sitance as log10
  gi$dist_log10 <- log10(InteractionSet::pairdist(gi, type = "mid"))

  # remove pairs with distance >= maxDist
  # this is only essential in case of non-zero length ranges in inGR
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

#' Add correlation of anchor signals to pairs of close genomic regions.
#'
#' This function adds a vector with correlation values for each input
#' interaction. Only works for input interaction within the given \code{maxDist}
#' distance.
#'
#' @param gi A sorted \code{\link[InteractionSet]{GInteractions}} object.
#' @param datcol a string matching an annotation column in \code{regions(gi)}.
#'   This collumn is assumed to hold the same number of values for each
#'   interaction as a \code{NumericList}.
#' @param colname A string that is used as columnname for the new column in
#'   \code{gi}.
#' @param maxDist maximal distance of pairs in bp as numeric. If maxDist=NULL,
#'   the maximal distance is computed from input interactions gi by
#'   \code{max(pairdist(gi))}.
#' @return A \code{\link[InteractionSet]{GInteractions}} similar to \code{gi}
#'   just wiht an additinoal column added.
#' @export
#' @import data.table
addCovCor <- function(gi, datcol, colname = "cor",
                           maxDist = NULL){

  # Algorithm to avoid comparisons of distal pairs
  # (0) define maxDist
  # (1) Group genome in overlapping bins of size 2*maxDist
  # (2) Run pairwise correlation for all ranges in each bin
  # (3) Combine correlations to data.frame with proper id1 and id2 in first
  #     columns
  # (4) Query fubak data.frame with input pairs
  #   /uses inner_join() from dplyr like in
  #  http://stackoverflow.com/questions/26596305/match-two-data-frames-based-on-multiple-columns

  # check input
  if ( any(is.na(GenomeInfoDb::seqlengths(gi))) ) {
    stop("gi object need seqlengths.")
  }
  if ( !methods::is(gi, "StrictGInteractions") ) {
    stop("gi should be of class StrictGInteractions")
  }

  #-----------------------------------------------------------------------------
  # (0) define maxDist
  #-----------------------------------------------------------------------------
  if (is.null(maxDist)) {

    # to avoid too many bins use minimal size of 10^6
    maxDist <- max(10^6, InteractionSet::pairdist(gi))

    if (maxDist < max(InteractionSet::pairdist(gi))) {
      stop(paste0("maxDist is smaller than maximal distance between",
                  "interactions in input gi."))
    }
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

  # get assign bins to regions
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

      m <- stats::cor(dat[,subIdx])

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
  matches <- corDT[gpDT, on = c("id1", "id2"), mult = "first"]

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
#' @param gi \code{\link{GInteractions}}.
#' @param scoreColname Character as name the metadata column in with motif score.
#' @return The same \code{\link{GInteractions}} as \code{gi} but with three
#'   additonal annotation columns.
#' @export
addMotifScore <- function(gi, scoreColname = "score"){

  anc1 <- InteractionSet::anchors(gi, type = "first", id = TRUE)
  anc2 <- InteractionSet::anchors(gi, type = "second", id = TRUE)

  ancScore <- S4Vectors::mcols(InteractionSet::regions(gi))[, scoreColname]

  S4Vectors::mcols(gi)[, "score_1"] <- ancScore[anc1]
  S4Vectors::mcols(gi)[, "score_2"] <- ancScore[anc2]
  S4Vectors::mcols(gi)[, "score_min"] <- apply(cbind(ancScore[anc1], ancScore[anc2]), 1, min)

  return(gi)
}

#' Prepares motif pairs and add genomic features.
#'
#' @param motifs \code{\link[GenomicRanges]{GRanges}} object with motif
#'   locations.
#' @inheritParams getCisPairs
#' @inheritParams addStrandCombination
#' @inheritParams addMotifScore
#'
#' @return An \code{\link[InteractionSet]{GInteractions}} object with moitf
#'   pairs and annotations of distance, strand orientation, and motif scores.
#'
#' @export
prepareCandidates <- function(motifs, maxDist = 10^6, scoreColname = "score"){

  # get pairs of motifs as GInteraction object
  gi <- getCisPairs(motifs, maxDist = maxDist)

  # add motif orientation combinations
  gi <- addStrandCombination(gi)

  # add motif score
  gi <- addMotifScore(gi, scoreColname = scoreColname)

  return(gi)
}

#' Add correlation of ChIP-seq coverage to motif pairs.
#'
#' This function first adds ChIP-seq signals along all regions of motif location
#' using the function \code{\link{addCovToGR}}. Than it calculates the
#' correlation of coverage for each input pair using the function \code{\link{addCovCor}}.
#' The pearson correlation value is added as new column metad data colum to the input interactions.
#'
#' @param gi \code{\link[InteractionSet]{GInteractions}} object.
#' @param bwFile File path or connection to BigWig file with ChIP-seq signals.
#' @param name Character indicating the sample name.
#' @inheritParams addCovToGR
#' @inheritParams addCovCor
#' @export
addCor <- function(gi, bwFile, name = "chip", window = 1000, binSize = 1){

  InteractionSet::regions(gi) <- addCovToGR(
    InteractionSet::regions(gi),
    bwFile,
    colname = name,
    window = window,
    binSize = binSize
    )

  # compute correlation of ChIP-seq profiles
  gi <- addCovCor(gi, datcol = name, colname = name)

  return(gi)
}

