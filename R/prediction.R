
#' Predict interaction probability using logistic regression model.
#'
#' @param data  A data.frame like object with predictor variables.
#' @param formula A \code{\link[stats]{formula}}. All predictor variables should
#'   be available in the data frame.
#' @param betas A vector with parameter estimates for predictor variables. They
#'   should be in the same order as variables in \code{formula}.
#'
#' @return A numeric vector with interaction probabilities for each observation
#'   in \code{df}. NAs are produced for NAs in \code{df}.
#'
predLogit <- function(data, formula, betas){

  # save option and set na.action to pass
  op <- options(na.action = "na.pass")
  on.exit(options(op), add = TRUE)

  # build model matrix
  mat <- stats::model.matrix(formula, data)

  pred <- as.numeric(betas %*% t(mat))
  pred <- boot::inv.logit(pred)

  return(pred)
}

#'Predict looping interactions.
#'
#'This function takes a \code{\link[InteractionSet]{GInteractions}} object with
#'candidate looping interactions. It should be annotated with features in
#'metadata columns. A logistic regression model is applied to predict looping
#'interaction probabilities.
#'
#'@param gi A \code{\link[InteractionSet]{GInteractions}} object with coverage
#'  correlation and genomic features in metadata columns. See
#'  \link{prepareCisPairs} and \link{addCor} to build it.
#'@param formula A \code{\link[stats]{formula}}. All predictor variables should
#'  be available in the in metadata columns of \code{gi}. If NULL, the following
#'  default formula is used: \code{~ dist + strandOrientation + score_min +
#'  chip}.
#'@param betas A vector with parameter estimates for predictor variables. They
#'  should be in the same order as variables in \code{formula}. Per default
#'  estimates of \code{\link{modelBest10Avg}} are used. See
#'  \code{?modelBest10Avg} for more detailed information on each parameter.
#'@param colname A \code{character} as column name of new metadata column in
#'  \code{gi} for predictions.
#'@param cutoff Numeric cutoff on prediction score. Only interactions with
#'  interaction probability >= \code{cutoff} are reported. If \code{NULL}, all
#'  input interactions are reported. Default is \code{\link{cutoffBest10}}, an
#'  optimal cutoff based on F1-score on 10 best performing transcription factor
#'  ChIP-seq data sets. See \code{?cutoffBest10} for more details.
#'
#'@return A \code{\link[InteractionSet]{GInteractions}} as \code{gi} with an
#'  additional metadata column holding the predicted looping probability.
#'
#'@seealso \code{\link{prepareCisPairs}}, \code{\link{addCor}}
#'
#' @examples
#'
#'# use example CTCF moitf location on human chromosome 22 with chip coverage
#'motifGR <- chromloop::motif.hg19.CTCF.chr22.cov
#'
#'# build candidate interactions
#'gi <- prepareCisPairs(motifGR, scoreColname = "sig")
#'
#'# add ChIP-seq signals correlation
#'gi <- addCovCor(gi)
#'
#'# predict chromatin looping interactions
#'loops <- predLoops(gi)
#'
#'# add prediction score for all candidates without filter
#'gi <- predLoops(gi, cutof = NULL)
#'
#'# add prediction score using custom column name
#'gi <- predLoops(gi, cutof = NULL, colname = "my_colname")
#'
#'# Filter loop predictions on custom cutoff
#'loops <- predLoops(gi, cutoff = 0.4)
#'
#'# predict chromatin looping interactions using custom model parameters
#'myParams <- c(-4, -5, -2, -1, -1, 5, 3)
#'loops <- predLoops(gi, betas = myParams)
#'
#'# predict chromatin loops using custom model formula and params
#'myFormula <- ~ dist + score_min
#'# define parameters for intercept, dist and motif_min
#'myParams <- c(-5, -4, 6)
#'loops <- predLoops(gi, formula = myFormula, betas = myParams)
#'
#'@import InteractionSet
#'@export
predLoops <- function(gi, formula = NULL, betas=NULL, colname = "pred",
                      cutoff = get("cutoffBest10")){
  # check arguments
  stopifnot(is(formula, "formula") | is.null(formula))
  stopifnot(is.numeric(betas) | is.null(betas))

  # if no formula is given, use default formula
  if (is.null(formula)) {
    formula <- stats::as.formula(
      " ~ dist + strandOrientation + score_min + cor_chip"
    )
  }

  # if no estimates are given, use the default parameters
  if (is.null(betas)) {

    # use default model
    defaultModel <- get("modelBest10Avg")

    betas <- defaultModel$estimate

  }

  # predict interaction
  pred <- predLogit(mcols(gi), formula, betas)

  # add prediction as new metadata column
  mcols(gi)[, colname] <- pred

  # test cutoff argument
  if (!is.null(cutoff)) {

    # filter interactions with by cutoff
    gi <- gi[!is.na(pred) & pred >= cutoff]
  }

  return(gi)
}
