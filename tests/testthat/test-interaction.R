context("interactions")


# Example input data ------------------------------------------------------

inGR <- GenomicRanges::GRanges(
  rep("chr1", 5),
  IRanges::IRanges(
    c(10, 20, 30, 100, 1000),
    c(15, 25, 35, 105, 1005)
  )
)

cov <- IRanges::RleList(list(
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

toySeqInfo <- GenomeInfoDb::Seqinfo(seqnames = c("chr1"),
                                        seqlengths = c(25),
                                        isCircular = c(FALSE),
                                        genome = "example")

toyCov <- IRanges::RleList(list(
  "chr1" = c(rep(0,4), 1:4, 4:1, rep(0,3), rep(1,4), 0, 1:3, 0, 0)
))

toyCovFile <- file.path(tempdir(), "toy.bw")
rtracklayer::export.bw(toyCov, toyCovFile)

toyGR <- GenomicRanges::GRanges(
  rep("chr1", 4),
  IRanges::IRanges(
    c(1, 5, 20, 14),
    c(4, 8, 23, 17)
  ),
  strand = c("+", "+", "+", "-"),
  score = c(5, 4, 6, 7),
  seqinfo = toySeqInfo
)

toyGI <- InteractionSet::GInteractions(
  c(1, 2, 2),
  c(4, 3, 4),
  toyGR,
  mode = "strict"
)


# Tests -------------------------------------------------------------------



test_that("getCisPairs works on small example dataset", {

  gi <- getCisPairs(inGR, 50)

  expect_equal(length(gi), 3)

  a1 <- InteractionSet::anchors(gi, type = "first", id = TRUE)
  a2 <- InteractionSet::anchors(gi, type = "second", id = TRUE)
  expect_true(all(a1 < a2))
  expect_true(all(InteractionSet::pairdist(gi) <= 50))

})

test_that("getCisPairs returns divergent pairs", {

  gr <- GenomicRanges::GRanges(
    rep("chr1", 4),
    IRanges::IRanges(
      c(5, 15, 1, 20),
      c(8, 18, 4, 23)
    ),
    strand = c("+", "+", "-", "-"),
    seqinfo = toySeqInfo
  )

  gi <- getCisPairs(gr, maxDist=10)

  gi <- addStrandCombination(gi)

  # test that idx 1 and 3 are contained
  anc1 <- InteractionSet::anchors(gi, type = "first", id = TRUE)
  anc2 <- InteractionSet::anchors(gi, type = "second", id = TRUE)

  expect_true(any(anc1 == 1 & anc2 == 3))

  int1_3 <- which(anc1 == 1 & anc2 == 3)

  expect_equal(gi$strandOrientation[int1_3], "divergent")

})

test_that("coverage is added to gi on small example dataset", {

  gi <- getCisPairs(inGR, 50)

  regGR <- InteractionSet::regions(gi)

  InteractionSet::regions(gi) <- addCovToGR(
    regGR,
    bwFile,
    window = 10,
    bin_size = 1,
    colname = "cov")

})

test_that("getCisPairs works with whole CTCF motif data", {

  skip("skipt test on whole CTCF motif data set for time.")

  # use internal motif data
  motifGR <- motif.hg19.CTCF

  # get all pairs within 1Mb
  gi <- getCisPairs(motifGR, 1e6)

  expect_equal(InteractionSet::regions(gi), motifGR)
  expect_true(max(InteractionSet::pairdist(gi)) <= 1e6)

})

test_that("applyToCloseGI runs on toy example dataset", {


  InteractionSet::regions(toyGI) <- addCovToGR(InteractionSet::regions(toyGI),
                                            toyCovFile,
                                            window = 4,
                                            bin_size = 1,
                                            colname = "cov")

  ncolBefore <- ncol(S4Vectors::mcols(toyGI))

  toyGI <- applyToCloseGI(toyGI, "cov", fun = cor, colname = "value")

  # check that exactly one column is added
  expect_equal(ncol(S4Vectors::mcols(toyGI)), ncolBefore + 1)
  expect_equal(names(S4Vectors::mcols(toyGI)), "value")

  # check that all cor values are correct

  covLst <- InteractionSet::regions(toyGI)$cov
  pairIdx <- cbind(
    InteractionSet::anchors(toyGI, id = TRUE, type = "first"),
    InteractionSet::anchors(toyGI, id = TRUE, type = "second")
  )
  expect_warning(
    corMat <- cor(t(as.matrix(covLst)))
  )
  manualCor <- corMat[pairIdx]

  expect_equal(toyGI$value, manualCor)

})

test_that("interactions can be annotated with Hi-C loops", {

  # skip("Skipt test on large files")

  # use internal motif data on chr22
  motifGR <- motif.hg19.CTCF[GenomeInfoDb::seqnames(motif.hg19.CTCF) == "chr22"]

  # get all pairs within 1Mb
  gi <- getCisPairs(motifGR, 1e6)

  # parse loops from internal data
  exampleLoopFile <- system.file(
    "extdata",
    "GM12878_HiCCUPS.chr22_1-18000000.loop.txt",
    package = "chromloop")

  loopGI <- parseLoopsRao(exampleLoopFile, seqinfo = GenomeInfoDb::seqinfo(gi))

  # hits <- InteractionSet::findOverlaps(gi, loopGI, type="within")

  ovl <- InteractionSet::countOverlaps(gi, loopGI, type = "within") >= 1

  S4Vectors::mcols(gi)[,"loop"] <- ovl

})

test_that("addInteractionSupport works with toy example", {

  toySupport <- InteractionSet::GInteractions(
    GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 4)),
    GenomicRanges::GRanges("chr1", IRanges::IRanges(15, 20))
  )

  toyGI <- addInteractionSupport(toyGI, toySupport)

  expect_equal(toyGI$loop == "Loop", c(TRUE, FALSE, FALSE))

})

test_that("addStrandCombination works for straned and unstraed ranges", {

  skip("TODO: implement test")

  toyGI <- addStrandCombination(toyGI)

})

test_that("addMotifScore works on toy example", {

  toyGI <- addMotifScore(toyGI, colname = "score")

  expect_true(all(c("score_1", "score_2", "score_min") %in% names(mcols(toyGI))))
  expect_equal(toyGI$score_1, toyGR$score[anchors(toyGI, "first", id = TRUE)])
  expect_equal(toyGI$score_min, apply(cbind(toyGI$score_1, toyGI$score_2), 1, min))

})

