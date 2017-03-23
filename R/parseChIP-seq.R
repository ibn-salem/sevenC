
#'  Parse a bigWig file as RleList coverage object
#'
#' The bigWig files contains read counts (or other dense, continouse data)
#' along the genome.
#' The bigWig format is described here:
#' \url{https://genome.ucsc.edu/goldenpath/help/bigWig.html}
#'
#' @param inFile Input file path or connection. See \code{con} paramter in
#' \code{\link[rtracklayer]{import}} function.
#' @param seqInfo A \code{\link{seqinfo}} object defining the reference genome.
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

  # # set seqinfo object
  # GenomeInfoDb::seqlevels(bwGR) <- GenomeInfoDb::seqlevels(seqInfo)
  # GenomeInfoDb::seqinfo(bwGR) <- seqInfo
  #
  # # compute coverage by using the score as weight
  # bwCov = GenomicRanges::coverage(bwGR, weight="score")

  # return(bwCov)
  return(cov)
}


getOutOfBound <- function(gr){

  # check if extended regions are out of chromsome space
  outOfBoundIdx <- GenomicRanges:::get_out_of_bound_index(gr)

  # save length of left and right out-of-bound lengths
  outGR <- gr[outOfBoundIdx]

  outDF <- data.frame(
    idx = outOfBoundIdx,
    left = ifelse(GenomicRanges::start(outGR) <= 0,
                  abs(GenomicRanges::start(outGR)) + 1,
                  0),
    right = ifelse(GenomicRanges::end(outGR) > GenomeInfoDb::seqlengths(outGR),
                   abs(GenomicRanges::end(outGR) - GenomeInfoDb::seqlengths(outGR)),
                   0)
  )

  return(outDF)

}


#' Parse coverage for specific regions from bigWig file.
#' @param inFile The connection from which data is loaded. If this is a
#'   character vector, it is assumed to be a filename and a corresponding file
#'   connection is created
#' @param ... Other parameters to pass to \code{\link[rtracklayer]{import.bw}}.
#' @value RleList or other object defined by \code{as}.
parseBigWigCov <- function(inFile, as = "RleList", ...){

  cov = rtracklayer::import.bw(con = inFile, as = as, ...)

  return(cov)
}

#' Add coverage to regions in \code{\link[GenomicRanges]{GRanges}} object.
#'
#' This function adds a coverage vector to each range in a genomic ranges
#' object. The coverage is reported for a fixed-sized window around the region
#' center and is reversed in case of negative strand of the region.
#'
#' @param gr \code{\link[GenomicRanges]{GRanges}} object with genomic regions
#' @param bwFile File path or connection to BigWig file with coverage to parrse
#'   from.
#' @param window the window size arund the center of ranges in \code{gr}.
#' @param bin_size size of bins to which the coverage values are combined. This
#'   is not implemented yet.
#' @return \code{\link[GenomicRanges]{GRanges}} as input but with an additional
#'   meta column containing the coverage values for each region.
#' @export
addCovToGR <- function(gr, bwFile, window=1000, bin_size=1, colname="cov"){

  # get windows around gr
  suppressWarnings(
    ancWin <- GenomicRanges::resize(gr, width = window, fix = "center")
    )

  outDF <- getOutOfBound(ancWin)

  # if (length(outOfBoundIdx) > 0) {
  #   stop("Windows around regions extend out of chromosomal bounds.")
  # }

  # # check taht window reg does not extend range of coverage data
  # covRange <- sapply(cov, length)
  # winRange <- range(range(ancWin))

  # trim ranges to fint within chromosomes
  ancWin <- GenomicRanges::trim(ancWin)
  
  message("INFO: Start reading coverage from fiel: ", bwFile, " ...")
  # get numeric with coverage of each region
  covGR <- rtracklayer::import.bw(
    bwFile,
    selection = rtracklayer::BigWigSelection(ancWin),
    as = "GRanges")
  message("INFO: Finished reading coverage from fiel: ", bwFile)
  
  # update covGR with seqinfo to allow subsetting with ancWin
  GenomeInfoDb::seqlevels(covGR) <- GenomeInfoDb::seqlevels(ancWin) 
  GenomeInfoDb::seqinfo(covGR) <- GenomeInfoDb::seqinfo(ancWin)
  
  covRle <- GenomicRanges::coverage(covGR, weight=covGR$score)
    
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

  # reverse coverage vector for regons on minus strand
  # TODO

  # add as additional column to GRanges object
  S4Vectors::mcols(gr)[,colname] <- covAnc

  return(gr)

}
