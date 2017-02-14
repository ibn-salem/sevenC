context("interactions")

test_that("getCisPairs works on small example dataset", {

  inGR <- GenomicRanges::GRanges(
    rep("chr1", 5),
    IRanges::IRanges(
      c(10, 20, 30, 100, 1000),
      c(15, 25, 35, 105, 1005)
    )
  )

  gi <- getCisPairs(inGR, 50)

  expect_equal(length(gi), 6)

  a1 <- InteractionSet::anchors(gi, type="first", id=TRUE)
  a2 <- InteractionSet::anchors(gi, type="second", id=TRUE)
  expect_true(all(a1 < a2))
  expect_true(all(InteractionSet::pairdist(gi) <= 50))

})
