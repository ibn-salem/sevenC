context("parseChIP-seq")

exampleBigWig <- system.file("extdata", "GM12878_Stat1.chr22_1-18000000.bigWig",
                             package = "chromloop")

toySeqInfo <- GenomeInfoDb::Seqinfo(seqnames = c("chr1", "chr22"),
                                    seqlengths = c(10^8, 10^8),
                                    isCircular = c(FALSE, FALSE),
                                    genome = "toy")

exampleGR <- GenomicRanges::GRanges(
  rep("chr22", 3),
  IRanges::IRanges(
    18000000 - c(1000, 500, 100),
    18000000 - c(991, 491, 91)
  ),
  seqinfo = GenomeInfoDb::Seqinfo(seqnames = "chr22",
                                  seqlengths = 10^8,
                                  isCircular = FALSE,
                                  genome = "toy_hr22")
)

test_path <- system.file("tests", package = "rtracklayer")
test_bw <- file.path(test_path, "test.bw")
# GRanges object with 9 ranges and 1 metadata column:
#   seqnames       ranges strand |     score
# <Rle>    <IRanges>  <Rle> | <numeric>
#   [1]     chr2 [   1,  300]      * |        -1
# [2]     chr2 [ 301,  600]      * |     -0.75
# [3]     chr2 [ 601,  900]      * |      -0.5
# [4]     chr2 [ 901, 1200]      * |     -0.25
# [5]     chr2 [1201, 1500]      * |         0
# [6]    chr19 [1501, 1800]      * |      0.25
# [7]    chr19 [1801, 2100]      * |       0.5
# [8]    chr19 [2101, 2400]      * |      0.75
# [9]    chr19 [2401, 2700]      * |         1
# -------
#   seqinfo: 2 sequences from an unspecified genome

# cov <- IRanges::RleList(list(
#   "chr1" = rep(c(0, 1, 0), each = 20),
#   "chr22" = rep(c(0, 5, 1, 5, 0), each = 10)
# ))

test_gr <- GenomicRanges::GRanges(
  rep("chr2", 3),
  IRanges::IRanges(
    c(350, 900, 1300),
    c(350, 900, 1300)
  ),
  seqinfo = GenomeInfoDb::Seqinfo(seqnames = c("chr2", "chr19"),
                                  seqlengths = c(3000, 3000),
                                  isCircular = c(FALSE, FALSE),
                                  genome = "example")
)


test_that("parseBigWigCov returns RleList of proper dimenstion", {

  cov <- parseBigWigCov(test_bw)

})

test_that("parseBigWigToRle on example dataset", {

  bwCov <- parseBigWigToRle(exampleBigWig, toySeqInfo)

  expect_true(length(bwCov) > 0)

})

test_that("parseBigWigToRle can parse only selected regions in example data", {

  bwCov <- parseBigWigToRle(exampleBigWig, toySeqInfo)

  expect_true(length(bwCov) > 0)

})


test_that("addCovToGR works on example dataset", {

  grCov <- addCovToGR(
    test_gr,
    test_bw, window = 10,
    bin_size = 1,
    colname = "covTest")

  expect_equal(length(test_gr), length(grCov))
  expect_true(all(all(!is.na(grCov$covTest))))
  expect_equal(length(S4Vectors::mcols(grCov)[, "covTest"]), length(test_gr))
  expect_equal(as.vector(sum(grCov$covTest[1])), 10 * -0.75)

})

test_that("addCovToGR handles out of chromosome coordinates", {

  out_gr <- GenomicRanges::GRanges(
    rep("chr2", 3),
    IRanges::IRanges(
      c(2, 700, 2999),
      c(2, 700, 2999)
    ),
    seqinfo = GenomeInfoDb::Seqinfo(seqnames = c("chr2", "chr19"),
                                    seqlengths = c(3000, 3000),
                                    isCircular = c(FALSE, FALSE),
                                    genome = "example")
  )

  grCov <- addCovToGR(out_gr, test_bw, window = 10, bin_size = 1,
                      colname = "covTest")

  expect_true(!is.null(grCov$covTest))
  expect_equal(grCov$covTest[[1]], c(rep(NA, 4), rep(-1, 6)))
  expect_true(all(sapply(grCov$covTest, length) == 10))
})

test_that("addCovToGR handles input ranges without seqinfo object", {

  out_gr <- GenomicRanges::GRanges(
    rep("chr2", 3),
    IRanges::IRanges(
      c(2, 700, 2999),
      c(2, 700, 2999)
    )
  )

  grCov <- addCovToGR(out_gr, test_bw, window = 10, bin_size = 1,
                      colname = "covTest")

  expect_true(!is.null(grCov$covTest))
  expect_equal(grCov$covTest[[1]], c(rep(NA, 4), rep(-1, 6)))
  expect_true(all(sapply(grCov$covTest, length) == 10))
  expect_true(!any(is.na(grCov$covTest[[3]])))

})


test_that("addCovToGR reversed coverage for regions on negative strand.", {

  stranedGR <- test_gr
  GenomicRanges::strand(stranedGR) <- c("+", "-", "-")

  grCov <- addCovToGR(test_gr, test_bw, window = 10, bin_size = 1,
                      colname = "covTest")
  stranedCov <- addCovToGR(stranedGR, test_bw, window = 10, bin_size = 1,
                           colname = "covTest")

  expect_equal(stranedCov$covTest[[2]], rev(grCov$covTest[[2]]) )
  expect_equal(stranedCov$covTest[[1]], grCov$covTest[[1]])
})

test_that("addCovToGR handles bin_size and window paramter correctly.", {

  grCov <- addCovToGR(
    test_gr,
    test_bw, window = 100,
    bin_size = 10,
    colname = "covTest")

  expect_equal(length(test_gr), length(grCov))
  expect_true(all(all(!is.na(grCov$covTest))))
  expect_equal(length(S4Vectors::mcols(grCov)[, "covTest"]), length(test_gr))
  expect_equal(length(grCov$covTest[[1]]), 10)
})

