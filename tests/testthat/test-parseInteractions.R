context("parseInteractions")


exampleLoopFile <- system.file("extdata", "GM12878_HiCCUPS.chr22_1-18000000.loop.txt", package="chromloop")

toySeqInfo <- GenomeInfoDb::Seqinfo(seqnames=c("chr1", "chr22"),
                                    seqlengths=c(10^8, 10^8),
                                    isCircular=c(FALSE, FALSE),
                                    genome="toy")

test_that("parseLoopsRao works with example file", {

  gi <- parseLoopsRao(exampleLoopFile)

  df <- read.table(exampleLoopFile, header=TRUE)

  expect_equal(length(gi),  nrow(df))
  expect_true(all(
    GenomicRanges::seqnames(InteractionSet::regions(gi)) == "chr22"
    ))
})

test_that("parseLoopsRao works toy genome as seqinfo", {

  gi <- parseLoopsRao(exampleLoopFile, seqinfo=toySeqInfo)

  expect_identical(GenomeInfoDb::seqinfo(gi), toySeqInfo)

})
