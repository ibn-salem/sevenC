context("parseChIP-seq")

exampleBigWig <- system.file("extdata", "GM12878_Stat1.chr22_1-18000000.bigWig",
                             package = "sevenC")

toySeqInfo <- Seqinfo(seqnames = c("chr1", "chr22"),
                                    seqlengths = c(10^8, 10^8),
                                    isCircular = c(FALSE, FALSE),
                                    genome = "toy")

exampleGR <- GRanges(
  rep("chr22", 3),
  IRanges(
    18000000 - c(1000, 500, 100),
    18000000 - c(991, 491, 91)
  ),
  seqinfo = Seqinfo(seqnames = "chr22",
                                  seqlengths = 10^8,
                                  isCircular = FALSE,
                                  genome = "toy_hr22")
)

test_path <- system.file("tests", package = "rtracklayer")
test_bw <- file.path(test_path, "test.bw")
test_wig <- file.path(test_path, "step.wig")

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

# cov <- RleList(list(
#   "chr1" = rep(c(0, 1, 0), each = 20),
#   "chr22" = rep(c(0, 5, 1, 5, 0), each = 10)
# ))

test_gr <- GRanges(
  rep("chr2", 3),
  IRanges(
    c(350, 900, 1300),
    c(350, 900, 1300)
  ),
  seqinfo = Seqinfo(seqnames = c("chr2", "chr19"),
                                  seqlengths = c(3000, 3000),
                                  isCircular = c(FALSE, FALSE),
                                  genome = "example")
)


test_that("addCovToGR works on example dataset and gives warning on windows", {

  if (.Platform$OS.type == 'windows') {

    expect_warning(
      addCovToGR(
        test_gr,
        test_bw, window = 10,
        binSize = 1,
        colname = "covTest")
    )

  } else {

    grCov <- addCovToGR(
      test_gr,
      test_bw, window = 10,
      binSize = 1,
      colname = "covTest")

    expect_equal(length(test_gr), length(grCov))
    expect_true(all(all(!is.na(grCov$covTest))))
    expect_equal(length(mcols(grCov)[, "covTest"]), length(test_gr))
    expect_equal(as.vector(sum(grCov$covTest[1])), 10 * -0.75)
  }
})

test_that("addCovToGR handles out of chromosome coordinates", {

  # bigWig parsing fails on windows in rtracklayer::import.bw()
  if (.Platform$OS.type != 'windows') {

    out_gr <- GRanges(
      rep("chr2", 3),
      IRanges(
        c(2, 700, 2999),
        c(2, 700, 2999)
      ),
      seqinfo = Seqinfo(seqnames = c("chr2", "chr19"),
                                      seqlengths = c(3000, 3000),
                                      isCircular = c(FALSE, FALSE),
                                      genome = "example")
    )

    grCov <- addCovToGR(out_gr, test_bw, window = 10, binSize = 1,
                        colname = "covTest")

    expect_true(!is.null(grCov$covTest))
    expect_equal(grCov$covTest[[1]], c(rep(NA, 4), rep(-1, 6)))
    expect_true(all(sapply(grCov$covTest, length) == 10))
  }
})

test_that("addCovToGR handles input ranges without seqinfo object", {

  skip("Currently not supported.")

  # bigWig parsing fails on windows in rtracklayer::import.bw()
  if (.Platform$OS.type != 'windows') {

    out_gr <- GRanges(
      rep("chr2", 3),
      IRanges(
        c(2, 700, 2999),
        c(2, 700, 2999)
      )
    )

    grCov <- addCovToGR(out_gr, test_bw, window = 10, binSize = 1,
                        colname = "covTest")

    expect_true(!is.null(grCov$covTest))
    expect_equal(grCov$covTest[[1]], c(rep(NA, 4), rep(-1, 6)))
    expect_true(all(sapply(grCov$covTest, length) == 10))
    expect_true(!any(is.na(grCov$covTest[[3]])))
  }
})

test_that("addCovToGR handles chromosomes not contained in bigWig file", {

  # bigWig parsing fails on windows in rtracklayer::import.bw()
  if (.Platform$OS.type != 'windows') {

    out_gr <- GRanges(
      c(rep("chr2", 3), "chrNew"),
      IRanges(
        c(2, 700, 2999, 10),
        c(2, 700, 2999, 20)
      ),
      seqinfo = Seqinfo(
        c("chr2", "chrNew"),
        c(5000, 100),
        c(FALSE, FALSE),
        genome = "chrNewGenome"
      )
    )

    # expect Warning that chrNew is not contained in bigWig file
    testthat::expect_warning(
      grCov <- addCovToGR(out_gr, test_bw, window = 10, binSize = 1,
                        colname = "covTest")
    )
    expect_true(!is.null(grCov$covTest))
    expect_equal(grCov$covTest[[1]], c(rep(NA, 4), rep(-1, 6)))
    expect_true(all(sapply(grCov$covTest, length) == 10))
    expect_true(!any(is.na(grCov$covTest[[3]])))
  }
})

test_that("addCovToGR reversed coverage for regions on negative strand.", {

  stranedGR <- test_gr
  strand(stranedGR) <- c("+", "-", "-")

  # bigWig parsing fails on windows in rtracklayer::import.bw()
  if (.Platform$OS.type != 'windows') {
    grCov <- addCovToGR(test_gr, test_bw, window = 10, binSize = 1,
                        colname = "covTest")
    stranedCov <- addCovToGR(stranedGR, test_bw, window = 10, binSize = 1,
                             colname = "covTest")

    expect_equal(stranedCov$covTest[[2]], rev(grCov$covTest[[2]]) )
    expect_equal(stranedCov$covTest[[1]], grCov$covTest[[1]])
  }

})

test_that("addCovToGR handles binSize and window paramter correctly.", {

  # bigWig parsing fails on windows in rtracklayer::import.bw()
  if (.Platform$OS.type != 'windows') {

    grCov <- addCovToGR(
      test_gr,
      test_bw, window = 100,
      binSize = 10,
      colname = "covTest")

    expect_equal(length(test_gr), length(grCov))
    expect_true(all(all(!is.na(grCov$covTest))))
    expect_equal(length(mcols(grCov)[, "covTest"]), length(test_gr))
    expect_equal(length(grCov$covTest[[1]]), 10)
  }
})


test_that("addCovToGR works with chr22 example data", {

  # use internal motif data
  motifGR <- sevenC::motif.hg19.CTCF.chr22

  # bigWig parsing fails on windows in rtracklayer::import.bw()
  if (.Platform$OS.type != 'windows') {

    motifGR <- addCovToGR(motifGR, exampleBigWig)

    expect_true("chip" %in% names(mcols(motifGR)))
  }
})

test_that("addCovToGR works with wig file", {


  # build custom gr with proper chroms as in step.wig
  chr19_gr <- GRanges(
    rep("chr19", 3),
    IRanges(
      c(59104901, 59105901, 59106301),
      c(59104911, 59105911, 59106311)
    ),
    seqinfo = Seqinfo(seqnames = c("chr18", "chr19"),
                      seqlengths = c(60000000, 60000000),
                      isCircular = c(FALSE, FALSE),
                      genome = "chr19_gr")
  )

  # run with wig file
  grCov <- addCovToGR(
    chr19_gr,
    test_wig, window = 10,
    binSize = 1,
    colname = "covTest")

    expect_equal(length(chr19_gr), length(grCov))
    expect_true(all(all(!is.na(grCov$covTest))))
    expect_equal(length(mcols(grCov)[, "covTest"]), length(chr19_gr))
    # compare to actual step.wig file from rtracklayer package
    expect_equal(grCov$covTest[[1]][1], 12.5)

})
