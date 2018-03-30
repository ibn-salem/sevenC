#' Build a \code{\link{GInteractions}} object with all pairs of input
#' \code{\link[GenomicRanges]{GRanges}} within a given distance.
#'
#' Distance is calculated from the center of input regions.
#'
#' @param inGR \code{\link[GenomicRanges]{GRanges}} object of genomic regions.
#'   The ranges should be sorted according to chr, strand, and start position.
#'   Use \code{\link[BiocGenerics]{sort}} to sort it.
#' @param maxDist maximal distance in base-pairs between pairs of ranges as
#'   single  numeric value.
#' @return A \code{\link[InteractionSet]{GInteractions}} object with all pairs
#'   within the given distance.
#' @examples
#'# build example GRanges as input
#'inGR <- GRanges(
#' rep("chr1", 5),
#' IRanges(
#'   c(10, 20, 30, 100, 1000),
#'   c(15, 25, 35, 105, 1005)
#' )
#')
#'
#'# get all pairs within 50 bp
#'gi <- getCisPairs(inGR, maxDist = 50)
#'
#'# getCisPiars returns a StrictGInteractions object
#'class(gi)
#'
#'# The input regions are accessibly via regions()
#'regions(gi)
#'
#' @import InteractionSet
#' @importFrom BiocGenerics sort
#' @importFrom GenomicRanges resize findOverlaps
#' @export
getCisPairs <- function(inGR, maxDist = 1e6){

  # check that input GRanges object is sorted
  if (!all(inGR == sort(inGR))) stop("Input ranges inGR need to be sorted. Use
                                     sort(inGR) to sort it")

  # get center postions of each input range
  posGR <- resize(inGR, width = 1, fix = "center")

  # calculate overlap all possible gene pairs within maxDist bp
  hits <- findOverlaps(posGR,
                      maxgap = maxDist,
                      drop.redundant = TRUE,
                      drop.self = TRUE,
                      ignore.strand = TRUE)

  # build IntractionSet object
  gi <- GInteractions(
    queryHits(hits),
    subjectHits(hits),
    inGR,
    mode = "strict"
  )

  # sort gi
  gi <- sort(gi)


  # add distance
  gi$dist <- pairdist(gi, type = "mid")

  # remove pairs with distance >= maxDist
  # this is only essential in case of non-zero length ranges in inGR
  gi <- gi[gi$dist <= maxDist]

  return(gi)
}


#' returns indices of columns with non-zero variance
#'
#' @param dat data.frame or matrix
#' @return column indices of columns with non-zero variance
noZeroVar <- function(dat) {
  out <- apply(dat, 2, function(x) length(unique(x)))
  which(out > 1)
}

#'Add correlation of anchor signals to pairs of close genomic regions.
#'
#'This function adds a vector with correlation values for each input
#'interaction. Only works for input interaction within the given \code{maxDist}
#'distance. Note, this function does not work on windows because reading of
#'bigWig files is currently not supported on windows.
#'
#'@param gi A sorted \code{\link[InteractionSet]{GInteractions}} object.
#'@param datacol a string matching an annotation column in \code{regions(gi)}.
#'  This column is assumed to hold the same number of values for each
#'  interaction as a \code{NumericList}.
#'@param colname A string that is used as columnname for the new column in
#'  \code{gi}.
#'@param maxDist maximal distance of pairs in bp as numeric. If maxDist=NULL,
#'  the maximal distance is computed from input interactions gi by
#'  \code{max(pairdist(gi))}.
#'@param use an optional character string giving a method for computing
#'  covariances in the presence of missing values. See \code{\link[stats]{cor}}
#'  for more details.
#'@param method a character string indicating which correlation coefficient (or
#'  covariance) is to be computed. One of "pearson" (default), "kendall", or
#'  "spearman": can be abbreviated. See \code{\link[stats]{cor}} for more
#'  details.
#'
#'@return A \code{\link[InteractionSet]{GInteractions}} similar to \code{gi}
#'  just with an additional column added.
#' @examples
#'if (.Platform$OS.type != "windows") {
#'
#'   # use internal motif data on chromosome 22
#'   motifGR <- sevenC::motif.hg19.CTCF.chr22
#'
#'   # use example bigWig file
#'  exampleBigWig <- system.file("extdata",
#'      "GM12878_Stat1.chr22_1-30000000.bigWig", package = "sevenC")
#'
#'   # add coverage from bigWig file
#'   motifGR <- addCovToGR(motifGR, exampleBigWig)
#'
#'   # get all pairs within 1Mb
#'   gi <- getCisPairs(motifGR, 1e5)
#'
#'   # compute correaltion of coverge for each pair
#'   gi <- addCovCor(gi)
#'
#'   # addCovCor adds a new metadata column:
#'   mcols(gi)
#'
#'}
#'@import data.table
#'@import InteractionSet
#'@importFrom BiocGenerics start
#'@importFrom GenomeInfoDb seqlengths seqinfo keepSeqlevels seqnames
#'@importFrom GenomicRanges GRanges findOverlaps slidingWindows
#'@importFrom S4Vectors mcols mcols<- queryHits subjectHits
#'@importFrom methods is
#'@export
addCovCor <- function(gi, datacol = "chip", colname = "cor_chip",
                      maxDist = NULL, use = "everything", method = "pearson"){



  # Algorithm to avoid comparisons of distal pairs
  # (0) define maxDist
  # (1) Group genome in overlapping bins of size 2*maxDist
  # (2) Run pairwise correlation for all ranges in each bin
  # (3) Combine correlations to data.frame with proper id1 and id2 in first
  #     columns
  # (4) Query fubak data.frame with input pairs

  # check input
  if ( any(is.na(seqlengths(gi))) ) {
    stop("gi object need seqlengths.")
  }
  if ( !is(gi, "StrictGInteractions") ) {
    stop("gi should be of class StrictGInteractions")
  }

  #-----------------------------------------------------------------------------
  # (0) define maxDist
  #-----------------------------------------------------------------------------
  if (is.null(maxDist)) {

    # to avoid too many bins use minimal size of 10^6
    maxDist <- max(10e6, pairdist(gi))

    if (maxDist < max(pairdist(gi))) {
      stop(paste0("maxDist is smaller than maximal distance between",
                  "interactions in input gi."))
    }
  }

  #-----------------------------------------------------------------------------
  # (1) group ranges in by bins
  #-----------------------------------------------------------------------------

  # create GRanges object for entire genome
  genomeGR <- GRanges(seqinfo(gi))

  # tile genoe in overlapping bins of with 2*maxDist
  binGR <- unlist(slidingWindows(genomeGR, 2 * maxDist, maxDist))

  # get assign bins to regions
  hits <- findOverlaps(binGR, regions(gi))

  #-----------------------------------------------------------------------------
  # (2) compute pairwise correlatin for all ranges in each bin
  #-----------------------------------------------------------------------------

  covList <- mcols(regions(gi))[, datacol]
  datamat <- as.matrix(covList)

  corMatList <- lapply(1:length(binGR), function(i){

    # get regions in this bin
    regIdx <- subjectHits(hits)[queryHits(hits) == i]

    if (length(regIdx) == 1) {
      dat <- cbind(datamat[regIdx, ])
    }else{
      dat <- t(datamat[regIdx, ])
    }

    # get indices with non-zero variance (they casue warning and NA in cor())
    subIdx <- noZeroVar(dat)

    n <- length(subIdx)

    # compute pairwise correlations for all regions in this bin
    if (n != 1) {

      m <- stats::cor(dat[, subIdx], use = use, method = method)

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

  corDT <- data.table::rbindlist(corMatList)

  #-----------------------------------------------------------------------------
  # (4) Query with input pairs
  #-----------------------------------------------------------------------------

  names(corDT) <- c("id1", "id2", "val")
  data.table::setkeyv(corDT, cols = c("id1", "id2"))

  # convert gp into data.table and set keys to id1 and id2 columns
  gpDT <- data.table::data.table(
    id1 = anchors(gi, type = "first", id = TRUE),
    id2 = anchors(gi, type = "second", id = TRUE),
    key = c("id1", "id2")
  )


  matches <- corDT[gpDT, on = c("id1", "id2"), mult = "first"]

  mcols(gi)[, colname] <- matches$val

  return(gi)
}


#' Add column to \code{\link{GInteractions}} with overlap support.
#'
#' See overlap methods in \code{\link{InteractionSet}} package for more details
#' on the overlap calculations: \code{?overlapsAny}
#'
#' @param gi \code{\link{GInteractions}} object
#' @param subject another \code{\link{GInteractions}} object
#' @param colname name of the new annotation column in \code{gi}.
#' @param ... additional arguments passed to \code{\link[IRanges]{overlapsAny}}.
#' @return \code{\link{InteractionSet}} \code{gi} as input but with additional
#'   annotation column \code{colname} indicating whether each interaction
#'   is supported by \code{subject} or not.
#' @examples
#'
#' # build example GRanges as anchors
#'anchorGR <- GRanges(
#'  rep("chr1", 4),
#'  IRanges(
#'    c(1, 5, 20, 14),
#'    c(4, 8, 23, 17)
#'  ),
#'  strand = c("+", "+", "+", "-"),
#'  score = c(5, 4, 6, 7)
#')
#'
#'
#'# build example GIntreaction object
#'gi <- GInteractions(
#'  c(1, 2, 2),
#'  c(4, 3, 4),
#'  anchorGR,
#'  mode = "strict"
#')
#'
#'# build exapple support GInteractions object
#'exampleSupport <- GInteractions(
#'     GRanges("chr1", IRanges(1, 4)),
#'     GRanges("chr1", IRanges(15, 20))
#')
#'
#'# add support
#'gi <- addInteractionSupport(gi, subject = exampleSupport)
#'
#'# Use colname argument to add support to differnt metadata column name
#'gi <- addInteractionSupport(gi, subject = exampleSupport, colname = "example")
#'
#' @import InteractionSet
#' @importFrom IRanges overlapsAny
#' @export
addInteractionSupport <- function(gi, subject, colname = "loop", ...){

  ol <- overlapsAny(gi, subject, ...)

  mcols(gi)[, colname] <- factor(ol,
                                          c(FALSE, TRUE),
                                          c("No loop", "Loop"))

  return(gi)
}


#' Add combination of anchor strand orientation.
#'
#' Each anchor region has a strand that is \code{'+'} or \code{'-'}. Therefore,
#' the each interaction between two regions has one of the following strand
#' combinations: "forward", "reverse", "convergent", or "divergent". Unstranded
#' ranges, indicated by \code{*}, are treated as positive strand.
#' @param gi \code{\link{GInteractions}}
#' @param colname name of the new column that is created in \code{gi}.
#'
#' @return The same \code{\link{GInteractions}} as \code{gi} but with an
#'   additional column indicating the four possible combinations of strands
#'   "forward", "reverse", "convergent", or "divergent".
#'
#' @examples
#'
#'# build example GRanges as anchors
#'anchorGR <- GRanges(
#'  rep("chr1", 4),
#'  IRanges(
#'    c(1, 5, 20, 14),
#'    c(4, 8, 23, 17)
#'  ),
#'  strand = c("+", "+", "+", "-"),
#'  score = c(5, 4, 6, 7)
#')
#'
#'
#'# build example GIntreaction object
#'gi <- GInteractions(
#'  c(1, 2, 2),
#'  c(4, 3, 4),
#'  anchorGR,
#'  mode = "strict"
#')
#'
#'# add combination of anchor strands as new metadata column
#'gi <- addStrandCombination(gi)
#'
#'# build small matrix to check strand combination
#'cbind(
#' as.character(strand(anchors(gi, "first"))),
#' as.character(strand(anchors(gi, "second"))),
#' mcols(gi)[, "strandOrientation"]
#')
#'
#' @import InteractionSet
#' @importFrom BiocGenerics strand
#' @export
addStrandCombination <- function(gi, colname = "strandOrientation"){

  anc1 <- anchors(gi, type = "first", id = TRUE)
  anc2 <- anchors(gi, type = "second", id = TRUE)

  # because interactions in GInteraction are sorted that + strand comes first
  # we need to make sure that the more upstream anchor region is taken as first
  ancStart <- start(regions(gi))
  anc1First <- ancStart[anc1] <= ancStart[anc2]

  first <- ifelse(anc1First, anc1, anc2)
  second <- ifelse(anc1First, anc2, anc1)

  ancStrand <- strand(regions(gi))

  sameStrand <- ancStrand[first] == ancStrand[second]
  firstPos <- as.character(ancStrand[first]) %in% c("+", "*")

  mcols(gi)[sameStrand & firstPos, colname] <- "forward"
  mcols(gi)[sameStrand & !firstPos, colname] <- "reverse"
  mcols(gi)[!sameStrand & firstPos, colname] <- "convergent"
  mcols(gi)[!sameStrand & !firstPos, colname] <- "divergent"

  return(gi)
}


#' Add motif score of anchors.
#'
#' If each anchor region (motif) has a score as annotation column, this function
#' adds two new columns named "score_1" and "score_2" with the scores of the
#' first and the second anchor region, respectively. Additionally, a column
#' named "score_min" is added with holds for each interaction the minimum of
#' "score_1" and "score_2".
#' @param gi \code{\link{GInteractions}}.
#' @param scoreColname Character as name the metadata column in with motif
#'   score.
#' @return The same \code{\link{GInteractions}} as \code{gi} but with three
#'   additional annotation columns.
#'
#' @examples
#'
#'# build example GRanges as anchors
#'anchorGR <- GRanges(
#'  rep("chr1", 4),
#'  IRanges(
#'    c(1, 5, 20, 14),
#'    c(4, 8, 23, 17)
#'  ),
#'  strand = c("+", "+", "+", "-"),
#'  score = c(5, 4, 6, 7)
#')
#'
#'
#'# build example GIntreaction object
#'gi <- GInteractions(
#'  c(1, 2, 2),
#'  c(4, 3, 4),
#'  anchorGR,
#'  mode = "strict"
#')
#'
#'# add add motif score
#'gi <- addMotifScore(gi, scoreColname = "score")
#'
#' @import InteractionSet
#' @export
addMotifScore <- function(gi, scoreColname = "score"){

  # check arguments:
  stopifnot(scoreColname %in% names(mcols(regions(gi))))

  anc1 <- anchors(gi, type = "first", id = TRUE)
  anc2 <- anchors(gi, type = "second", id = TRUE)

  ancScore <- mcols(regions(gi))[, scoreColname]

  mcols(gi)[, "score_1"] <- ancScore[anc1]
  mcols(gi)[, "score_2"] <- ancScore[anc2]
  mcols(gi)[, "score_min"] <- apply(
    cbind(ancScore[anc1], ancScore[anc2]), 1, min
    )

  return(gi)
}

#' Prepares motif pairs as \code{\link[InteractionSet]{GInteractions}} and add
#' genomic features.
#'
#' @param motifs \code{\link[GenomicRanges]{GRanges}} object with motif
#'   locations.
#' @inheritParams getCisPairs
#' @inheritParams addStrandCombination
#' @inheritParams addMotifScore
#'
#' @return An \code{\link[InteractionSet]{GInteractions}} object with motif
#'   pairs and annotations of distance, strand orientation, and motif scores.
#'
#' @examples
#'
#'# build example GRanges as anchors
#'anchorGR <- GRanges(
#'  rep("chr1", 4),
#'  IRanges(
#'    c(1, 5, 20, 14),
#'    c(4, 8, 23, 17)
#'  ),
#'  strand = c("+", "+", "+", "-"),
#'  score = c(5, 4, 6, 7)
#')
#'
#'
#'
#'# prepare candidates
#'gi <- prepareCisPairs(anchorGR)
#'
#'
#'# prepare candidates using a mimial distance of 10 bp
#'gi <- prepareCisPairs(anchorGR, maxDist = 10)
#'
#'# prepare candidates using an alternative score value in anchors
#'anchorGR$myScore <- rnorm(length(anchorGR))
#'gi <- prepareCisPairs(anchorGR, scoreColname = "myScore")
#'
#' @import InteractionSet
#' @export
prepareCisPairs <- function(motifs, maxDist = 1e6, scoreColname = "score"){

  # get pairs of motifs as GInteraction object
  gi <- getCisPairs(motifs, maxDist = maxDist)

  # add motif orientation combinations
  gi <- addStrandCombination(gi)

  # add motif score
  gi <- addMotifScore(gi, scoreColname = scoreColname)

  return(gi)
}

#'Add correlation of ChIP-seq coverage to motif pairs.
#'
#'This function first adds ChIP-seq signals along all regions of motif location
#'using the function \code{\link{addCovToGR}}. Than it calculates the
#'correlation of coverage for each input pair using the function
#'\code{\link{addCovCor}}. The Pearson correlation coefficient is added as new
#'metadata column to the input interactions. Note, this function does not work
#'on windows because reading of bigWig files is currently not supported on
#'windows.
#'
#'@param gi \code{\link[InteractionSet]{GInteractions}} object.
#'@param bwFile File path or connection to BigWig file with ChIP-seq signals.
#'@param name Character indicating the sample name.
#'@inheritParams addCovToGR
#'@inheritParams addCovCor
#'
#'@return An \code{\link[InteractionSet]{GInteractions}} object like \code{gi}
#'  with a new metadata column \code{colname} holding Pearson correlation
#'  coefficient of ChIP-seq signals for each anchor pair.
#'
#'@examples
#'if (.Platform$OS.type != "windows") {
#'
#'  # use example bigWig file of ChIP-seq signals on human chromosome 22
#'  exampleBigWig <- system.file("extdata",
#'  "GM12878_Stat1.chr22_1-30000000.bigWig", package = "sevenC")
#'
#'  # use example CTCF moitf location on human chromosome 22
#'  motifGR <- sevenC::motif.hg19.CTCF.chr22
#'
#'  # build candidate interactions
#'  gi <- prepareCisPairs(motifGR)
#'
#'
#'  # add ChIP-seq signals correlation
#'  gi <- addCor(gi, exampleBigWig)
#'
#'  # use an alternative metadata column name for ChIP-seq correlation
#'  gi <- addCor(gi, exampleBigWig, name = "Stat1")
#'
#'
#'  # add ChIP-seq correlation for signals signals in windows of 500bp around
#'  # motif centers
#'  gi <- addCor(gi, exampleBigWig, window = 500)
#'
#'
#'  # add ChIP-seq correlation for signals in bins of 10 bp
#'  gi <- addCor(gi, exampleBigWig, binSize = 10)
#'
#'}
#'@import InteractionSet
#'@export
addCor <- function(gi, bwFile, name = "chip", window = 1000, binSize = 1){

  regions(gi) <- addCovToGR(
    regions(gi),
    bwFile,
    colname = name,
    window = window,
    binSize = binSize
    )

  # compute correlation of ChIP-seq profiles
  gi <- addCovCor(gi, datacol = name, colname = paste0("cor_", name))

  return(gi)
}
