context("predictions")

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

