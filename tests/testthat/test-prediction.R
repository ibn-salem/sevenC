context("predictions")


# prepare gi from chr22 exampple data set
exampleGI <- prepareCisPairs(motif.hg19.CTCF.chr22.cov, scoreColname = "sig")

exampleGI <- addCovCor(exampleGI)

# Tests -------------------------------------------------------------------

test_that("predLogit gives same result as predict.glm()", {


  formula <- am ~ wt + cyl + hp

  mod <- glm(formula, mtcars, family = binomial())
  betas <- coef(mod)

  pred <- predLogit(mtcars, betas, formula = formula)
  pred_glm <- predict.glm(mod, newdata = mtcars, type = "response")
  names(pred_glm) <- NULL

  expect_equal(pred, pred_glm)
})


test_that("predLogit works on chr22 example data", {

  # run prediction
  pred <- predLogit(
    mcols(exampleGI),
    stats::as.formula(" ~ dist + strandOrientation + score_min + cor_chip"),
    modelBest10Avg$estimate
  )

  expect_true(is.numeric(pred))
  expect_equal(length(pred), length(exampleGI))
  expect_true(all(pred <= 1 & pred >= 0, na.rm = TRUE))

})

test_that("predLoops works with default parameters and without cutoff", {

  # run prediction
  predGI <- predLoops(exampleGI, cutoff = NULL)

  expect_true(is.numeric(predGI$pred))
  expect_equal(length(predGI),  length(exampleGI))
  expect_true(all(predGI$pred <= 1 & predGI$pred >= 0, na.rm = TRUE))
  expect_equal(
    mcols(predGI)[, -ncol(mcols(predGI))],
    mcols(exampleGI))
})


test_that("predLoops works with default cutoff", {

  # run prediction
  predGI <- predLoops(exampleGI)

  expect_true(is.numeric(predGI$pred))
  expect_true(length(predGI) <= length(exampleGI))
  expect_true(all(predGI$pred >= cutoffBest10))
  expect_equal(
    ncol(mcols(predGI)),
    ncol(mcols(exampleGI)) + 1)
})
