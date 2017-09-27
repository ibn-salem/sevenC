
#'  Parse a bigWig file as RleList coverage object
#'
#' The bigWig files contains read counts (or other dense, continouse data)
#' along the genome.
#' The bigWig format is described here:
#' \url{https://genome.ucsc.edu/goldenpath/help/bigWig.html}
#'
#' @param inFile Input file path or connection. See \code{con} paramter in
#' \code{\link[rtracklayer]{import}} function.
#' @param seqInfo A \code{\link[GenomeInfoDb]{seqinfo}} object defining the reference genome.
#' @param selectionGR \code{\link[GenomicRanges]{GRanges}} object with regions for
#'  which the coverage is selected.
#' @return An \code{\link[IRanges]{RleList}} object with density values for each
#' position in the genome.
#' @export
parseBigWigToRle <- function(inFile, seqInfo, selectionGR = NULL){

  # if selection is given, use it.
  if (is.null(selectionGR)) {

    cov = rtracklayer::import.bw(con = inFile, as = "RleList")

  } else {

    # bwSelection <- rtracklayer::BigWigSelection(selectionGR)

    # parse file as GRange object
    cov = rtracklayer::import.bw(con = inFile,
                                  which = selectionGR,
                                  # selection = bwSelection,
                                  as = "RleList")
  }

  return(cov)
}


#'Parse BAM file as RleList coverage object.
#'
#'See http://bioconductor.org/packages/release/bioc/html/Rsamtools.html
#'
#'@param inFile a BAM file
#'@param seqInfo a \code{\link[GenomeInfoDb]{seqinfo}} object.
#'@param selectionGR \code{\link[GenomicRanges]{GRanges}} object with regions for
#'  which the coverage is selected.
#'@return \code{\link[S4Vectors]{Rle-class}} coverage object. See
#'  \code{\link[IRanges]{coverage}}.
#'
#'@importFrom Rsamtools ScanBamParam scanBamHeader
parseBAMToRle <- function(inFile, seqInfo, selectionGR = NULL){

  # inFile <- system.file("extdata", "ex1.bam", package = "Rsamtools")
  inFile <- "../chromloop_analysis/data/ENCODE/bam/wgEncodeSydhTfbsGm12878Rad21IggrabAlnRep1.bam"
  # selectionGR <- GRanges(paste("seq", 1:2, sep=""), IRanges(1, 1e7))
  selectionGR <- motif.hg19.CTCF

  seqInfo <- seqinfo(selectionGR)
  seqlevelsStyle(selectionGR) <- "UCSC"
  seqlevels(selectionGR)

  param <- ScanBamParam(
    what = c( "pos" , "qwidth" ),
    which = selectionGR,
    flag = scanBamFlag(isUnmappedQuery = FALSE)
    )
  h <- scanBamHeader(inFile)

  # x <- scanBam(inFile, param = param, seqinfo = seqinfo(selectionGR))[[1]]

  coverage(IRanges(x[["pos"]], width=x[["qwidth"]]))
}

#' Get out of chromosomal bound ranges.
#'
#' @param gr A \code{GRanges} object.
#' @return A \code{data.frame} with rows for each range in \code{gr} that
#'   extends out of chromosomes. The first colum holds the index of the range in
#'   \code{gr}, the second the size of the overlap to the left of the chromsome
#'   and the third the sieze of the overlap to the right of the chromosme.
getOutOfBound <- function(gr){

  # check if extended regions are out of chromsome space
  outOfBoundIdx <- GenomicRanges:::get_out_of_bound_index(gr)
  negIdx <- which( GenomicRanges::start(gr) <= 0 )

  outIdx <- union(outOfBoundIdx, negIdx)

  # save length of left and right out-of-bound lengths
  outGR <- gr[outIdx]

  chromEnds <- GenomeInfoDb::seqlengths(outGR)[
    as.vector(GenomeInfoDb::seqnames(outGR))]

  outDF <- data.frame(
    idx = outIdx,
    left = ifelse(
      GenomicRanges::start(outGR) <= 0,
      abs(GenomicRanges::start(outGR)) + 1,
      0),
    right = ifelse(
      !is.na(chromEnds) &
      GenomicRanges::end(outGR) >
        chromEnds,
      abs(GenomicRanges::end(outGR) - chromEnds),
      0)
  )

  return(outDF)

}


#' Parse coverage for specific regions from bigWig file.
#' @param inFile The connection from which data is loaded. If this is a
#'   character vector, it is assumed to be a filename and a corresponding file
#'   connection is created
#' @param ... Other parameters to pass to \code{\link[rtracklayer]{import.bw}}.
#' @param as Specifies the class of the return object. See \code{as} argument of
#'   \code{\link[rtracklayer]{import.bw}}.
#' @return RleList or other object defined by \code{as}.
parseBigWigCov <- function(inFile, as = "RleList", ...){

  cov = rtracklayer::import.bw(con = inFile, as = as, ...)

  return(cov)
}

#' Sliding mean over x of intervals of size k
#'
#' Source:
#' \url{http://stats.stackexchange.com/questions/3051/mean-of-a-sliding-window-in-r}
#'
#' @param x numeric vector
#' @param k interval size
#' @return numeric vector of length \code{length(x) / k}.
slideMean <- function(x, k){

  n <- length(x)
  spots <- seq(from = 1, to = n, by = k)
  result <- vector(length = length(spots))

  # for (i in 1:length(spots)) {
  #   result[i] <- mean(x[spots[i]:(spots[i] + k - 1)], na.rm = TRUE)
  # }

  result <- purrr::map_dbl(seq_along(spots), function(i){
    mean(x[seq(spots[i], (spots[i] + k - 1))], na.rm = TRUE)
  })

  return(result)

}


#' Add coverage to regions in \code{\link[GenomicRanges]{GRanges}} object.
#'
#' This function adds a vector of coverage (or any other signal in the input
#' bigWig file) to each range in a genomic ranges object. The coverage is
#' reported for a fixed-sized window around the region center. For regions with
#' negative strand, the coverage vector is reversed.
#'
#' @param gr \code{\link[GenomicRanges]{GRanges}} object with genomic regions
#' @param bwFile File path or connection to BigWig file with coverage to parrse
#'   from.
#' @param window the window size arund the center of ranges in \code{gr}.
#' @param bin_size size of bins to which the coverage values are combined.
#' @param colname name of the new colum that is created in \code{gr}.
#'
#' @return \code{\link[GenomicRanges]{GRanges}} as input but with an additional
#'   meta column containing the coverage values for each region.
#' @export
addCovToGR <- function(gr, bwFile, window=1000, bin_size=1, colname="cov"){

  # get windows around gr without warnding if ranges extend chromosome borders
  suppressWarnings(
    ancWin <- GenomicRanges::resize(gr, width = window, fix = "center")
  )

  # Because quering coverage result in error for ranges outisde chromosome we
  # need to get intervals defined outside of chromosomes
  outDF <- getOutOfBound(ancWin)

  # trim ranges to fint within chromosomes
  ancWin <- GenomicRanges::trim(ancWin)

  # trim start in case there is no seqinof object
  GenomicRanges::start(ancWin) <- ifelse(
    GenomicRanges::start(ancWin) > 0,
    GenomicRanges::start(ancWin),
    1)

  # get numeric with coverage of each region

  # define query region and trim seqinfo to avoid warning in rtracklayer::import.bw
  selectWin <- ancWin
  GenomeInfoDb::seqlevels(selectWin) <- as.character(
    unique(GenomeInfoDb::seqnames(selectWin))
    )
  GenomeInfoDb::seqinfo(selectWin) <- GenomeInfoDb::Seqinfo(
    GenomeInfoDb::seqlevels(selectWin)
    )
  selection <- rtracklayer::BigWigSelection(selectWin)

  message("INFO: Start reading coverage from file: ", bwFile, " ...")

  covGR <- rtracklayer::import.bw(
    bwFile,
    selection = selection,
    as = "GRanges",
    seqinfo = seqinfo(ancWin))
  message("INFO: Finished reading coverage from fiel: ", bwFile)

  # update covGR with seqinfo to allow subsetting with ancWin
  if ( !any(is.na(GenomeInfoDb::seqlengths(ancWin))) ) {
    GenomeInfoDb::seqlevels(covGR) <- GenomeInfoDb::seqlevels(ancWin)
    GenomeInfoDb::seqinfo(covGR) <- GenomeInfoDb::seqinfo(ancWin)
  }

  covRle <- GenomicRanges::coverage(covGR, weight = covGR$score)

  # get coverage as for all regions
  covAncRle <- covRle[ancWin]
  covAnc <- IRanges::NumericList(covAncRle)

  # add NAs for out of bound regions
  covAnc[outDF$idx] <- lapply(seq_along(outDF$idx), function(i){
    covAnc[outDF$idx][[i]] <- c(
      rep(NA, outDF[i, "left"]),
      covAnc[outDF$idx][[i]],
      rep(NA, outDF[i, "right"])
    )
  })

  # combine coverage for bins
  if (bin_size > 1) {
    covAnc <- IRanges::NumericList(lapply(covAnc, slideMean, k = bin_size))
  }


  # reverse coverage vector for regons on minus strand
  negStrand <- which(as.logical(GenomicRanges::strand(gr) == "-"))
  covAnc[negStrand] <- lapply(covAnc[negStrand], rev)


  # add as additional column to GRanges object
  S4Vectors::mcols(gr)[,colname] <- covAnc

  return(gr)

}


#' Parse signals along the genome for input regions
#'
#' @param files Character vector of bigWig files to be parsed
#' @param gr A \code{\link[GenomicRanges]{GRanges}} object with regions for which
#' the signals should be parased.
#' @param colData A data.frame like object for metadata of \code{files}.
#'
#' @return A \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}} object
# parseRangeSignals <- function(files, gr, colData){
#
#   # example code
#   # cd <- tibble(name = c("a", "b"))
#   cd <- tibble(name = c("a"))
#   nl <- NumericList(1:10, 1:10)
#   nlMat <- as.matrix(as.list(nl))
#
#   SummarizedExperiment(cbind(a=as.list(1:10), b=as.list(1:10)), colData = cd)
#   SummarizedExperiment(nlMat, colData = cd)
# }
