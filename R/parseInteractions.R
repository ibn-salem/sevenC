

#' Parse chromatin loops from Rao et al 2014 as an \code{\link{GInteractions}}
#' object.
#'
#' @param inFile input file with loops
#' @param ... additional arguments, that will be passed to
#'   \code{\link{GRanges()}} functions.
#' @return An \code{\link{GInteractions}} with loops from input file.
#' @export
parseLoopsRao <- function(inFile, ...){

  # parse IMR90 domains from Rao et al 2014:
  raoDF = utils::read.delim(inFile)

  # get chromsome column
  chr = paste0("chr", raoDF$chr1)

  # substract 1 bp from down coordinate to have inclusive interval ranges
  upAnchor = GenomicRanges::GRanges(
      chr,
      IRanges::IRanges(raoDF[,"x1"], raoDF[,"x2"] - 1),
      ...)

  downAnchor = GenomicRanges::GRanges(
    chr,
    IRanges::IRanges(raoDF[,"y1"], raoDF[,"y2"] - 1),
    ...)

  # build GInteractions
  gi <- InteractionSet::GInteractions(upAnchor, downAnchor, mode = "strict")

  # define and add additional annotations to each interaction
  annotationCols <- setdiff(names(raoDF),
                            c("chr1", "x1", "x2", "chr2", "y1", "y2"))

  S4Vectors::mcols(gi) <- raoDF[, annotationCols]

  return(gi)
}


#'Parse chromatin interactions from Tang et al 2015 as \code{GInteractions}.
#'
#'Reads pairwise ChIA-PET interaction from an input file.
#'
#'It reads files with the follwoing tab-delimited format:
#'
#'\tabular{lllllll}{ chr12\tab 48160351\tab 48161634\tab  chr12\tab 48230665\tab
#'48232848\tab 27 \cr chr7\tab 77284664\tab 77285815\tab chr7\tab  77388242\tab
#'77388928\tab 7 \cr chr4\tab 128459961\tab 128460166\tab chr4\tab
#'128508304\tab 128509082\tab 4 \cr }
#'
#'This file format was used for ChIA-PET interaction data by Tang et al. 2015
#'\url{http://dx.doi.org/10.1016/j.cell.2015.11.024}. The last column of input
#'file is added as annotation column with colname "score".
#'
#'@param inFile input file with loops
#'@param ... additional arguments, that will be passed to
#'  \code{\link{GRanges()}} functions.
#'@return An \code{\link{GInteractions}} with loops from input file.
#'@export
parseLoopsTang2015 <- function(inFile, ...){



  message("INFO: Parse ChiA-PET interactions from file: ", inFile)

  # parse IMR90 domains from Rao et al 2014:
  inDF = read.delim(inFile, header = FALSE)

  # substract 1 bp from down coordinate to have inclusive interval ranges
  upAnchor = GenomicRanges::GRanges(
    inDF[,1],
    IRanges::IRanges(inDF[,2], inDF[,3]), ...)

  downAnchor = GenomicRanges::GRanges(
    inDF[,4],
    IRanges::IRanges(inDF[,5], inDF[,6]), ...)


  # build GInteractions
  gi <- InteractionSet::GInteractions(upAnchor, downAnchor, mode = "strict")

  # annotate GInteractions object with score
  gi$score <- inDF[,7]

  return(gi)

}


