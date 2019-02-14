#'
#'@title Randomized Lasso 
#'
#'@description This function perform randomized lasso using the \code{glmnet} package. The 
#'function present in the \code{stabs} package that runs the lasso version was adapted for
#'the randomized lasso here. Randmized lasso stability selection uses this function repeatedly
#'to select predictors. 
#'
#'
#'@param x the predictor matrix 
#'@param y the response vector
#'@param q the average number fo selected variables
#'@param weakness the parameter essential for the randomized lasso. 
#'@param type 
#'
#'@return the regression output
#'
#'@author Dania Machlab
#'@export
glmnet.randomized_lasso <- function (x, y, q, weakness=1, type = c("conservative", "anticonservative"), ...) {
  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("Package ", sQuote("glmnet"), " needed but not available")
  if (is.data.frame(x)) {
    message("Note: ", sQuote("x"), " is coerced to a model matrix without intercept")
    x <- model.matrix(~. - 1, x)
  }
  if ("lambda" %in% names(list(...)))
    stop("It is not permitted to specify the penalty parameter ",
         sQuote("lambda"), " for lasso when used with stability selection.")
  type <- match.arg(type)
  # modify the function here to make it a randomized-lasso
  if (type == "conservative")
    fit <- suppressWarnings(glmnet::glmnet(x, y, pmax = q, penalty.factor = 1/runif(ncol(x), weakness,  1),  ...))
  if (type == "anticonservative")
    fit <- glmnet::glmnet(x, y, dfmax = q - 1, penalty.factor = 1/runif(ncol(x), weakness,  1), ...)
  selected <- predict(fit, type = "nonzero")
  selected <- selected[[length(selected)]]
  ret <- logical(ncol(x))
  ret[selected] <- TRUE
  names(ret) <- colnames(x)
  cf <- fit$beta
  sequence <- as.matrix(cf != 0)
  return(list(selected = ret, path = sequence))
}



#'@title Randomized Lasso Stability Selection
#'
#'@description This function runs randomized lasso stability selection as presented by Meinshausen (2010) 
#'and Buehlmann and with the improved error bounds introduced by Shah and Samworth (2013). The function 
#'uses the \code{stabsel} function from the \code{stabs} package, but using the randomized lasso version.
#'
#'@param x the predictor matrix 
#'@param y the response vector
#'@param weakness the parameter essential for the randomized lasso. 
#'@param cutoff value between 0 and 1 which is the cutoff for the selection probability. Any variable
#'with a selection probability that is higher than the set cutoff will be selected.
#'@param PFER the absolute number of false positives that we allow for in the final list of 
#'selected variables. This comes from the Meinshausen and Buehlmann paper.  
#'
#'@return a \code{stabsel} object
#'
#'@details Randomized lasso stability selection runs randomized lasso several times on several subsamples of 
#'the response variable and predictor matrix (stability selection). N/2 elements from the response variable are
#'randomly chosen, where N is the length of the vector. Their corresponsing section of the predictor matrix is 
#'also chosen, and the \code{glmnet.randomized_lasso} function is applied. This is done multiple times, and 
#'results in selection probabilities for each predictor. The probability of a specific predictor is the number of 
#'times it was selected divided by the total number of subsamples that were made (total number of times the 
#'regression was performed). 
#'
#'We make use of the \code{stabs} package that implements lasso stability selection, and adapt it to run the 
#'randomized lasso stability selection. The output is an object of type \code{stabsel}. 
#'
#'@author Dania Machlab
#'@export
randomized_stabsel <- function(x=x, y=y, weakness=0.8, cutoff=0.8, PFER=2, ...) {
  stabs::stabsel(x=x, y=y, fitfun=glmnet.randomized_lasso, args.fitfun=list(weakness=weakness), cutoff=cutoff, PFER=PFER, ...)
}





