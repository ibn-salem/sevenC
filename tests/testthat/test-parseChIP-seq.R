context("parseChIP-seq")


test_that("parseBigWigToRle on example dataset", {

  exampleBigWig <- system.file("extdata", "GM12878_Stat1.chr22_1-18000000.bigWig", package="chromloop")

  toySeqInfo <- GenomeInfoDb::Seqinfo(seqnames=c("chr1", "chr22"),
                                      seqlengths=c(10^8, 10^8),
                                      isCircular=c(FALSE, FALSE),
                                      genome="toy")

  bwCov <- parseBigWigToRle(exampleBigWig, toySeqInfo)

  expect_true(length(bwCov) > 0)
  expect_equal(length(bwCov), length(toySeqInfo))
  expect_equal(length(as.vector(bwCov$chr22)), 100000000)

})
