
#' Parse chromatin loops from Rao et al 2014 as
#' \code{\link[InteractionSet]{StrictGInteractions}}.
#'
#' @param inFile input file with loops
#' @param ... additional arguments, that will be passed to
#'   \code{\link[GenomicRanges]{GRanges}} functions.
#' @return \code{\link{GInteractions}} with loops from input file.
#' @export
parseLoopsRao <- function(inFile, ...){

  # parse input file
  raoDF = as.data.frame(readr::read_tsv(inFile))

  # get chromsome column
  chr = paste0("chr", raoDF$chr1)

  # substract 1 bp from end coordinate to have inclusive interval ranges
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
#'  \code{\link[GenomicRanges]{GRanges}} functions.
#'@return An \code{\link{GInteractions}} with loops from input file.
#'@export
parseLoopsTang2015 <- function(inFile, ...){

  message("INFO: Parse ChiA-PET interactions from file: ", inFile)

  # parse IMR90 domains from Rao et al 2014:
  inDF = readr::read_tsv(inFile, col_names = FALSE)

  # create ranges by adding +1 to start coordintate to convert from 0-based to
  # to 1-based coordinates.
  upAnchor = GenomicRanges::GRanges(
    inDF[[1]],
    IRanges::IRanges(inDF[[2]] + 1, inDF[[3]]), ...)

  downAnchor = GenomicRanges::GRanges(
    inDF[[4]],
    IRanges::IRanges(inDF[[5]] + 1, inDF[[6]]), ...)

  # build GInteractions
  gi <- InteractionSet::GInteractions(upAnchor, downAnchor, mode = "strict")

  # annotate GInteractions object with score
  gi$score <- inDF[[7]]

  return(gi)

}

#'Parse Capture Hi-C intreactions from Mifsud et al. 2015
#'
#' Capture Hi-C interactios from the study Mifsud et al. 2015 can be downloaded
#' from here
#' \url{http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2323/E-MTAB-2323.additional.1.zip}.
#' Each file in the .zip archive can be parsed with this function.
#'
#'@param inFile input file with interactions
#'@param ... additional arguments, that will be passed to
#'  \code{\link[GenomicRanges]{GRanges}} functions.
#'@return An \code{\link{GInteractions}} with interactions from input file.
#'@export
parseCaptureHiC <- function(inFile, ...){

  message("INFO: Parse interactions from file: ", inFile)

  inDF <- readr::read_tsv(inFile, col_names = TRUE)

  upAnchor = GenomicRanges::GRanges(
    inDF[[1]],
    IRanges::IRanges(inDF[[2]], inDF[[3]]),
    # Symbol = inDF[[4]],
    # Ensembl_Gene_ID = inDF[[5]],
    # expresssion_quartile = inDF[[6]],
    ...)

  downAnchor = GenomicRanges::GRanges(
    inDF[[7]],
    IRanges::IRanges(inDF[[8]], inDF[[9]]), ...)

  # only promoter-promoter files have the follwoing information
  # if ( ncol(inDF) > 11) {
  #   downAnchor$Symbol <- inDF[[10]]
  #   downAnchor$Ensembl_Gene_ID <- inDF[[11]]
  #   downAnchor$expresssion_quartile <- inDF[[12]]
  # }

  # build GInteractions
  gi <- InteractionSet::GInteractions(upAnchor, downAnchor)

  # annotate GInteractions object with score
  gi$raw_count <- inDF$`raw count`
  gi$log_observed_expected <- inDF$`log(observed/expected)`

  # turne into StrictGInteraction object
  gi <- methods::as(gi, "StrictGInteractions")

  return(gi)
}

