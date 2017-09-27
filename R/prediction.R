

#' Predict interaction probability using logistic regresstion model.
#'
#' @param data  A data.frame with predictor variables
#' @param formula A modelling formula. All predictor variables should be
#'   availabel in the data frame.
#' @param betas A vector with parameter estimates for predictor variables
#'
#' @return A numeric vector with interaction probabilites for each ovservation
#'   in \code{df}. NAs are produced for NAs in \code{df}.
#' @export
pred_logit <- function(data, formula, betas){

  # save option and set na.action to pass
  op <- options(na.action = 'na.pass')
  on.exit(options(op), add = TRUE)

  # build model matrix
  mat <- modelr::model_matrix(data, formula)

  pred <- as.numeric(betas %*% t(mat))
  pred <- boot::inv.logit(pred)

  return(pred)
}
