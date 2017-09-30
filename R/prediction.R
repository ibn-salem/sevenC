
#' Predict interaction probability using logistic regresstion model.
#'
#' @param data  A data.frame like object with predictor variables.
#' @param formula A \code{\link[stats]{formula}}. All predictor variables should
#'   be availabel in the data frame.
#' @param betas A vector with parameter estimates for predictor variables. They
#'   should be in the same order as variables in \code{formula}.
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

#' Predict looping interaction probability.
#'
#' @param gi A \code{\link[InteractionSet]{GInteractions}} object with coverage
#'   correlation and genomic features in metadata columns. See
#'   \link{prepareCandidates} and \link{addCor} to build it.
#' @param formula A \code{\link[stats]{formula}}. All predictor variables should
#'   be availabel in the in metadat colums of \code{gi}. If NULL, the following
#'   default formular is used: \code{~ dist + strandOrientation + score_min +
#'   chip}.
#' @param betas A vector with parameter estimates for predictor variables. They
#'   should be in the same order as variables in \code{formula}.
#' @param colname A \code{character} as column name of new metadata colum in
#'   \code{gi} for predictions.
#' @param cutoff Numeric cutoff on prediction score. Only interactions with
#'   interaction probability >= \code{cutoff} are reported. If \code{NULL}, all
#'   input interactions are reported. Default is \code{\link{cutoffBest10}}, an
#'   optimal cutoff based on F1-score on 10 best performing TF ChIP-seq data
#'   sets. See \code{?'cutoffBest10'} for more details.
#'
#' @return A \code{\link[InteractionSet]{GInteractions}} as \code{gi} with an
#'   additional metadata colum holdin the predicted looping probability.
#'
#' @seealso \code{\link{prepareCandidates}}, \code{\link{addCor}},
#'   \code{\link{pred_logit}}
#' @export
#'
predLoops <- function(gi, formula = NULL, betas = NULL, colname = "pred",
                      cutoff = cutoffBest10){

  # if no estimates are given, use the default parameters
  if (is.null(betas)) {
    betas <- modelBest10Avg$estimate
  }

  # if no formula is given, use default formula
  if (is.null(formula)) {
    formula <- stats::as.formula(
      " ~ dist + strandOrientation + score_min + chip"
    )
  }

  # predict interaction
  pred <- pred_logit(S4Vectors::mcols(gi), formula, betas)

  # add prediction as new metadata column
  S4Vectors::mcols(gi)[, colname] <- pred

  # test cutoff argument
  if (!is.null(cutoff)) {

    # filter interactions with by cutoff
    gi <- gi[!is.na(pred) & pred >= cutoff]
  }

  return(gi)
}
