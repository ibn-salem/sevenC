context("interactions")


# Example input data ------------------------------------------------------

inGR <- GRanges(
  rep("chr1", 5),
  IRanges(
    c(10, 20, 30, 100, 1000),
    c(15, 25, 35, 105, 1005)
  )
)

cov <- RleList(list(
  "chr1" = c(rep(c(0, 1, 0), each = 20), rep(0, 1000)),
  "chr22" = c(rep(c(0, 5, 1, 5, 0), each = 10), rep(0, 1000))
))

bwFile <- file.path(tempdir(), "test.bw")
rtracklayer::export.bw(cov, bwFile)


# Toy example data set ----------------------------------------------------


#'      4-|       XX
#'      3-|      XXXX            X
#'      2-|     XXXXXX          XX
#'      1-|    XXXXXXXX   XXXX XXX
#' chr1   |1...5....1....1....2....2
#'                  0    5    0    5
#'         0000123443210001111012300
#'         <--><-->     <-->  <-->
#' grWin     1   2        4     3
#' strand    +   +        -     +
#' score     5   4        6     7
#'
#' pairs:
#' 1 4
#' 2 4
#' 2 3

toySeqInfo <- Seqinfo(seqnames = c("chr1"),
                                        seqlengths = c(25),
                                        isCircular = c(FALSE),
                                        genome = "example")

toyCov <- RleList(list(
  "chr1" = c(rep(0, 4), 1:4, 4:1, rep(0, 3), rep(1, 4), 0, 1:3, 0, 0)
))

toyCovFile <- file.path(tempdir(), "toy.bw")
rtracklayer::export.bw(toyCov, toyCovFile)

toyGR <- GRanges(
  rep("chr1", 4),
  IRanges(
    c(1, 5, 20, 14),
    c(4, 8, 23, 17)
  ),
  strand = c("+", "+", "+", "-"),
  score = c(5, 4, 6, 7),
  seqinfo = toySeqInfo
)

toyGI <- GInteractions(
  c(1, 2, 2),
  c(4, 3, 4),
  toyGR,
  mode = "strict"
)


# Tests -------------------------------------------------------------------



test_that("getCisPairs works on small example dataset", {

  gi <- getCisPairs(inGR, 50)

  expect_equal(length(gi), 3)

  a1 <- anchors(gi, type = "first", id = TRUE)
  a2 <- anchors(gi, type = "second", id = TRUE)
  expect_true(all(a1 < a2))
  expect_true(all(pairdist(gi) <= 50))

})

test_that("getCisPairs returns divergent pairs", {

  gr <- GRanges(
    rep("chr1", 4),
    IRanges(
      c(5, 15, 1, 20),
      c(8, 18, 4, 23)
    ),
    strand = c("+", "+", "-", "-"),
    seqinfo = toySeqInfo
  )

  gi <- getCisPairs(gr, maxDist = 10)

  gi <- addStrandCombination(gi)

  # test that idx 1 and 3 are contained
  anc1 <- anchors(gi, type = "first", id = TRUE)
  anc2 <- anchors(gi, type = "second", id = TRUE)

  expect_true(any(anc1 == 1 & anc2 == 3))

  int13 <- which(anc1 == 1 & anc2 == 3)

  expect_equal(gi$strandOrientation[int13], "divergent")

})

test_that("coverage is added to gi on small example dataset", {

  gi <- getCisPairs(inGR, 50)

  regGR <- regions(gi)

  regions(gi) <- addCovToGR(
    regGR,
    bwFile,
    window = 10,
    binSize = 1,
    colname = "cov")

})

test_that("addCovCor works with chr22 example data", {

  exampleBigWig <- system.file("extdata",
                               "GM12878_Stat1.chr22_1-18000000.bigWig",
                               package = "chromloop")

  # use internal motif data
  motifGR <- chromloop::motif.hg19.CTCF.chr22

  motifGR <- addCovToGR(motifGR, exampleBigWig)

  # get all pairs within 1Mb
  gi <- getCisPairs(motifGR, 1e5)

  gi <- addCovCor(gi, datacol = "cov")

})


test_that("getCisPairs works with whole CTCF motif data", {

  # use internal motif data
  motifGR <- chromloop::motif.hg19.CTCF

  # get all pairs within 1Mb
  gi <- getCisPairs(motifGR, 1e6)

  expect_equal(regions(gi), motifGR)
  expect_true(max(pairdist(gi)) <= 1e6)

})

test_that("addCovCor runs on toy example dataset", {


  regions(toyGI) <- addCovToGR(regions(toyGI),
                                            toyCovFile,
                                            window = 4,
                                            binSize = 1,
                                            colname = "cov")

  ncolBefore <- ncol(mcols(toyGI))

  toyGI <- addCovCor(toyGI, "cov", colname = "value")

  # check that exactly one column is added
  expect_equal(ncol(mcols(toyGI)), ncolBefore + 1)
  expect_equal(names(mcols(toyGI)), "value")

  # check that all cor values are correct

  covLst <- regions(toyGI)$cov
  pairIdx <- cbind(
    anchors(toyGI, id = TRUE, type = "first"),
    anchors(toyGI, id = TRUE, type = "second")
  )
  expect_warning(
    corMat <- stats::cor(t(as.matrix(covLst)))
  )
  manualCor <- corMat[pairIdx]

  expect_equal(toyGI$value, manualCor)

})

test_that("interactions can be annotated with Hi-C loops", {

  # use internal motif data on chr22
  motifGR <- chromloop::motif.hg19.CTCF.chr22

  # get all pairs within 1Mb
  gi <- getCisPairs(motifGR, 1e6)

  # parse loops from internal data
  exampleLoopFile <- system.file(
    "extdata",
    "GM12878_HiCCUPS.chr22_1-18000000.loop.txt",
    package = "chromloop")

  loopGI <- parseLoopsRao(exampleLoopFile, seqinfo = seqinfo(gi))

  ovl <- countOverlaps(gi, loopGI, type = "within") >= 1

  mcols(gi)[,"loop"] <- ovl

})

test_that("addInteractionSupport works with toy example", {

  toySupport <- GInteractions(
    GRanges("chr1", IRanges(1, 4)),
    GRanges("chr1", IRanges(15, 20))
  )

  toyGI <- addInteractionSupport(toyGI, toySupport)

  expect_equal(toyGI$loop == "Loop", c(TRUE, FALSE, FALSE))

})

test_that("addStrandCombination works for straned and unstraed ranges", {

  toyGI <- addStrandCombination(toyGI)

  unstrandedGR <- regions(toyGI)
  strand(unstrandedGR) <- "*"

  unstrandedGI <- GInteractions(
    c(1, 2, 2),
    c(4, 3, 4),
    unstrandedGR,
    mode = "strict"
  )

  unstrandedGI <- addStrandCombination(unstrandedGI)

  expect_true("strandOrientation" %in% names(mcols(toyGI)))
  expect_true("strandOrientation" %in% names(mcols(unstrandedGI)))
  expect_true(any(toyGI$strandOrientation != unstrandedGI$strandOrientation))

})

test_that("addMotifScore works on toy example", {

  toyGI <- addMotifScore(toyGI, scoreColname = "score")

  expect_true(all(c("score_1", "score_2", "score_min") %in% names(mcols(toyGI))))
  expect_equal(toyGI$score_1, toyGR$score[anchors(toyGI, "first", id = TRUE)])
  expect_equal(toyGI$score_min, apply(cbind(toyGI$score_1, toyGI$score_2), 1, min))

})

test_that("prepareCandidates works on toy example data", {

  gi <- prepareCandidates(toyGR, 10)

  expect_equal(length(gi), 3)
  expect_true("strandOrientation" %in% names(mcols(gi)))
  expect_true("score_min" %in% names(mcols(gi)))

})

test_that("addCor works on toy example data", {

  gi <- addCor(toyGI, toyCovFile, name = "toy", window = 4, binSize = 1)

  expect_true("toy" %in% names(mcols(regions(gi))))
  expect_equal(ncol(mcols(gi)), ncol(mcols(toyGI)) + 1)
  expect_true("toy" %in% names(mcols(gi)))

})

