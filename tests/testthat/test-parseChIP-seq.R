context("parseChIP-seq")


exampleBigWig <- system.file("extdata", "GM12878_Stat1.chr22_1-18000000.bigWig", package="chromloop")

toySeqInfo <- GenomeInfoDb::Seqinfo(seqnames=c("chr1", "chr22"),
                                    seqlengths=c(10^8, 10^8),
                                    isCircular=c(FALSE, FALSE),
                                    genome="toy")


test_that("parseBigWigToRle on example dataset", {

  bwCov <- parseBigWigToRle(exampleBigWig, toySeqInfo)

  expect_true(length(bwCov) > 0)
  expect_equal(length(bwCov), length(toySeqInfo))
  expect_equal(length(as.vector(bwCov$chr22)), 100000000)

})


test_that("addCovToGR on example dataset", {

  cov <- IRanges::RleList(list(
    "chr1"=rep(c(0, 1, 0), each=20),
    "chr22"=rep(c(0, 5, 1, 5, 0), each=10)
  ))

  gr <- GenomicRanges::GRanges(
    rep("chr1", 3),
    IRanges::IRanges(
      c(10, 20, 30),
      c(10, 30, 35)
    ),
    seqinfo=GenomeInfoDb::Seqinfo(seqnames=c("chr1", "chr2"),
                                  seqlengths=c(200, 200),
                                  isCircular=c(FALSE, FALSE),
                                  genome="example")
  )

  grCov <- addCovToGR(gr, cov, window=10, bin_size=1, colname="covTest")

  expect_equal(length(gr), length(grCov))
  expect_true(all(all(!is.na(grCov$covTest))))
  expect_equal(length(S4Vectors::mcols(grCov)[, "covTest"]), length(gr))

  expect_equal(as.vector(sum(grCov$covTest[1])), 0)

})

