
#' Parse chromatin loops from Rao et al. 2014 as strict
#' \code{\link[InteractionSet:InteractionSet-class]{GInteractions}}.
#'
#' @param inFile input file with loops
#' @param ... additional arguments, that will be passed to
#'   \code{\link[GenomicRanges:GRanges-class]{GRanges}} functions.
#' @return \code{\link{GInteractions}} with loops from input file.
#'
#' @examples
#'
#' # use example loop file
#'exampleLoopFile <- system.file("extdata",
#'   "GM12878_HiCCUPS.chr22_1-30000000.loop.txt", package = "sevenC")
#'
#'# read loops form example file:
#'gi <- parseLoopsRao(exampleLoopFile)
#'
#'
#'# read loops with custom seqinfo object:
#'customSeqInfo <- Seqinfo(seqnames = c("chr1", "chr22"),
#'    seqlengths = c(10^8, 10^8), isCircular = c(FALSE, FALSE),
#'    genome = "custom")
#'gi <- parseLoopsRao(exampleLoopFile, seqinfo = customSeqInfo)
#'
#' @import InteractionSet
#' @importFrom IRanges IRanges
#' @export
parseLoopsRao <- function(inFile, ...){

  # parse input file
  raoDF <- as.data.frame(readr::read_tsv(
    inFile,
    col_types = readr::cols(
      chr1 = readr::col_character(),
      chr2 = readr::col_character())
    ))

  # get chromosome column
  chr_1 <- paste0("chr", as.character(raoDF$chr1))
  chr_2 <- paste0("chr", as.character(raoDF$chr2))

  # substract 1 bp from end coordinate to have inclusive interval ranges
  upAnchor <- GRanges(
      chr_1,
      IRanges(raoDF[, "x1"], raoDF[, "x2"] - 1),
      ...)

  downAnchor <- GRanges(
    chr_2,
    IRanges(raoDF[, "y1"], raoDF[, "y2"] - 1),
    ...)

  # build GInteractions
  gi <- GInteractions(upAnchor, downAnchor, mode = "strict")

  # define and add additional annotations to each interaction
  annotationCols <- setdiff(names(raoDF),
                            c("chr1", "x1", "x2", "chr2", "y1", "y2"))

  mcols(gi) <- raoDF[, annotationCols]

  return(gi)
}


#'Parse chromatin interactions from Tang et al. 2015 as \code{GInteractions}.
#'
#'Reads pairwise ChIA-PET interaction from an input file.
#'
#'It reads files with the following tab-delimited format:
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
#' @param inFile input file with loops
#' @param ... additional arguments, that will be passed to
#'  \code{\link[GenomicRanges:GRanges-class]{GRanges}} functions.
#' @return An \code{\link{GInteractions}} with loops from input file.
#'
#' @examples
#'exampleLoopTang2015File <- system.file("extdata",
#'    "ChIA-PET_GM12878_Tang2015.chr22_1-30000000.clusters.txt",
#'    package = "sevenC")
#'
#'gi <- parseLoopsTang(exampleLoopTang2015File)
#'
#'# read loops with custom seqinfo object:
#'customSeqInfo <- Seqinfo(seqnames = c("chr1", "chr22"),
#'    seqlengths = c(10^8, 10^8), isCircular = c(FALSE, FALSE),
#'    genome = "custom")
#'gi <- parseLoopsTang(exampleLoopTang2015File, seqinfo = customSeqInfo)
#'
#'
#' @import InteractionSet
#' @importFrom IRanges IRanges
#' @export
parseLoopsTang <- function(inFile, ...){

  message("INFO: Parse ChiA-PET interactions from file: ", inFile)

  # parse IMR90 domains from Rao et al. 2014:
  inDF <- readr::read_tsv(inFile, col_names = FALSE)

  # create ranges by adding +1 to start coordintate to convert from 0-based to
  # to 1-based coordinates.
  upAnchor <- GRanges(
    inDF[[1]],
    IRanges(inDF[[2]] + 1, inDF[[3]]), ...)

  downAnchor <- GRanges(
    inDF[[4]],
    IRanges(inDF[[5]] + 1, inDF[[6]]), ...)

  # build GInteractions
  gi <- GInteractions(upAnchor, downAnchor, mode = "strict")

  # annotate GInteractions object with score
  gi$score <- inDF[[7]]

  return(gi)

}

