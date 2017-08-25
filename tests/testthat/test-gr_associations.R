context("gr-associations")

# Toy example data set ----------------------------------------------------


#' chr1   |1...5....1....1....2....2....3
#'                  0    5    0    5....0
#' peaks     ==     ==
#' promoter   -              -
#'
#' loops: 1   |-|
#'        2         |--------|
#'        3      |-------------|
#'
#' anc        1 23  4        5 6

toySeqInfo <- GenomeInfoDb::Seqinfo(seqnames = c("chr1"),
                                        seqlengths = c(25),
                                        isCircular = c(FALSE),
                                        genome = "example")
toyPeaks <- GenomicRanges::GRanges(
  rep("chr1", 2),
  IRanges::IRanges(
    c(3, 10),
    c(4, 11)
  ),
  seqinfo = toySeqInfo
)


toyPromoter <- GenomicRanges::GRanges(
  rep("chr1", 2),
  IRanges::IRanges(
    c(4, 19),
    c(4, 19)
  ),
  seqinfo = toySeqInfo
)

toyAncGR <- GenomicRanges::GRanges(
  rep("chr1", 6),
  IRanges::IRanges(
    c(4, 7, 8, 11, 19, 21),
    c(4, 7, 8, 11, 19, 21)
  ),
  seqinfo = toySeqInfo
)

toyGI <- InteractionSet::GInteractions(
  c(1, 4, 3),
  c(2, 5, 6),
  toyAncGR,
  mode = "strict"
)

# Tests -------------------------------------------------------------------

test_that("linkRegions works on toy example dataset", {

  hits <- linkRegions(toyPeaks, toyPromoter, toyGI)

  # expected hits based on toy example data
  expHits <- data.frame(
    gr1 = c(1, 2),
    gr2 = c(1, 2),
    gi = c(NA, 2)
  )

  expect_equal(hits, expHits)
})

