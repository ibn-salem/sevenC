context("predictions")


# prepare gi from chr22 exampple data set
exampleGI <- prepareCandidates(motif.hg19.CTCF.chr22, scoreColname = "sig")

exampleBigWig <- system.file(
  "extdata", "GM12878_Stat1.chr22_1-18000000.bigWig", package = "chromloop")

exampleGI <- addCor(exampleGI, exampleBigWig)



# Tests -------------------------------------------------------------------

test_that("pred_logit gives same result as predict.glm()", {


  formula <- am ~ wt + cyl + hp

  mod <- glm(formula, mtcars, family = binomial())
  betas <- coef(mod)

  pred <- pred_logit(mtcars, betas, formula = formula)
  pred_glm <- predict.glm(mod, newdata = mtcars, type = "response")
  names(pred_glm) <- NULL

  expect_equal(pred, pred_glm)
})


test_that("pred_logit works on chr22 example data", {

  # run prediction
  pred <- pred_logit(
    S4Vectors::mcols(exampleGI),
    stats::as.formula(" ~ dist + strandOrientation + score_min + chip"),
    modelBest10Avg$estimate
  )

  expect_true(is.numeric(pred))
  expect_equal(length(pred), length(exampleGI))
  expect_true(all(pred <= 1 & pred >= 0, na.rm = TRUE))

})

test_that("predLoops works with default parameters", {

  # run prediction
  predGI <- predLoops(exampleGI)

  expect_true(is.numeric(predGI$pred))
  expect_equal(length(predGI), length(exampleGI))
  expect_equal(
      S4Vectors::mcols(predGI)[, -ncol(S4Vectors::mcols(predGI))],
      S4Vectors::mcols(exampleGI))
  expect_true(all(predGI$pred <= 1 & predGI$pred >= 0, na.rm = TRUE))

})

