context("gr-associations")

# Toy example data set ----------------------------------------------------


#' chr1   |1...5....1....1....2....2....3
#'                  0    5    0    5    0
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
    c(4, 6, 7, 10, 19, 21),
    c(4, 6, 7, 10, 19, 21)
  ),
  seqinfo = toySeqInfo
)

toyGI <- InteractionSet::GInteractions(
  c(1, 4, 3),
  c(2, 5, 6),
  toyAncGR,
  mode = "strict"
)

# strandedGI <- InteractionSet::GInteractions(
#   GenomicRanges::GRanges(
#     c("chr1", "chr1"),
#     IRanges::IRanges(
#       c(10, 40),
#       c(20, 50)),
#     strand = c("+", "+")),
#   GenomicRanges::GRanges(
#     c("chr1", "chr1"),
#     IRanges::IRanges(
#       c(15, 10),
#       c(25, 20)),
#     strand = c("+", "-")),
#     mode = "strict"
#   )

# Tests -------------------------------------------------------------------

test_that("linkRegions works on toy example dataset", {

  hits <- linkRegions(toyPeaks, toyPromoter, toyGI)

  # expected hits based on toy example data
  expHits <- data.frame(
    gr1 = c(1, 2),
    gr2 = c(1, 2),
    gi = c(NA, 2)
  )

  hitsInner3 <-  linkRegions(toyPeaks, toyPromoter, toyGI, inner_maxgap = 3)
  # expected hits based on toy example data
  expHitsInner3 <- data.frame(
    gr1 = c(1, 2, 2),
    gr2 = c(1, 2, 2),
    gi = c(NA, 2, 3)
  )

  expect_equal(hits, expHits)
  expect_equal(hitsInner3, expHitsInner3)
})

test_that("linkRegionsInLoops works on toy example dataset", {

  hits <- linkRegionsInLoops(toyPeaks, toyPromoter, toyGI)

  # expected hits based on toy example data
  expHits <- data.frame(
    gr1 = c(1, 1, 2, 2),
    gr2 = c(1, 1, 2, 2),
    gi = c(NA, 1, 2, 3)
  )

  expect_equal(hits, expHits)
})

# test_that("linkRegionsInLoops works on stranded anchros", {
#
#
#
# })


