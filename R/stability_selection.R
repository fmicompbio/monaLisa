#' @importFrom stabs stabsel
#' @importFrom glmnet glmnet predict.glmnet
#' @importFrom TFBSTools PWMatrixList
#' @importFrom stats model.matrix runif
NULL

#'@title Randomized Lasso
#'
#'@description This function performs randomized lasso using the \code{glmnet} package. The
#'function present in the \code{stabs} package that runs the lasso version was adapted for
#'the randomized lasso here. Randmized lasso stability selection uses this function repeatedly
#'to select predictors.
#'
#'
#'@param x the predictor matrix.
#'@param y the response vector.
#'@param q the average number of selected variables.
#'@param weakness parameter used randomized lasso (see details).
#'@param type parameter from \code{lars.lasso} function in \code{stabs}. It is a character vector specifying
#'how much the PFER should be controlled. If type is "conservative" then the number of selected variables per
#'subsample is <= q. If type is "anticonservative" then the number of selected variables per subsample is >= q.
#'By default the conservative type is chosen.
#'@param ... additional parameters for \code{glmnet}.
#'
#'@return the regression output which consists of a list of length 2. The list contains the following:
#'\itemize{
#'  \item selected - a logical vector of length equal to the total number of predictors. The predictors that were chosen have a value of TRUE.
#'  \item path - a logical matrix containing the regularization steps as columns and the predictors as rows. An entry of TRUE indicates selection.
#'}
#'
#'
#'@details This function is identical to \code{glmnet.lasso} from the \code{stabs} package. The only addition is the weakness parameter
#'which has been added when calling the \code{glmnet} function by setting penalty.factor = 1/runif(ncol(x), weakness, 1)
#'where ncol(x) is the number of predictors.
#'
#'@seealso \code{\link[stabs]{glmnet.lasso}} and \code{\link[glmnet]{glmnet}}
#'
#'@export
glmnet.randomized_lasso <- function(x, y, q, weakness=1, type = c("conservative", "anticonservative"), ...) {
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
  selected <- glmnet::predict.glmnet(fit, type = "nonzero")
  selected <- selected[[length(selected)]]
  ret <- logical(ncol(x))
  ret[selected] <- TRUE
  names(ret) <- colnames(x)
  cf <- fit$beta
  sequence <- as.matrix(cf != 0)
  return(list(selected = ret, path = sequence))
}



#' @title Randomized Lasso Stability Selection
#'
#' @description This function runs randomized lasso stability selection as presented by Meinshausen
#'   and Bühlmann (2010) and with the improved error bounds introduced by Shah and Samworth (2013). The function
#'   uses the \code{\link[stabs]{stabsel}} function from the \code{stabs} package, but implements the randomized lasso version.
#'
#' @param x the predictor matrix.
#' @param y the response vector.
#' @param weakness value between 0 and 1 (default = 0.8). 
#'   It affects how strict the method will be in selecting predictors. The closer it is to 0, 
#'   the more stringent the selection. A weakness value of 1 is identical to performing 
#'   lasso stability selection (not the randomized version).
#' @param cutoff value between 0 and 1 (default = 0.8) which is the cutoff for the selection probability. 
#'   Any variable with a selection probability that is higher than the set cutoff will be selected.
#' @param PFER integer (default = 2) representing the absolute number of false positives that we 
#'   allow for in the final list of selected variables. For details see Meinshausen and Bühlmann (2010).
#' @param ... additional parameters that can be passed on to \code{\link[glmnet]{glmnet}}.
#'
#' @details Randomized lasso stability selection runs a randomized lasso regression 
#'   several times on subsamples of the response variable and predictor matrix. 
#'   N/2 elements from the response variable are randomly chosen in each regression, 
#'   where N is the length of the vector. The corresponsing section of the predictor matrix is
#'   also chosen, and the \code{\link[monaLisa]{glmnet.randomized_lasso}} function is applied. 
#'   Stability selection results in selection probabilities for each predictor. 
#'   The probability of a specific predictor is the number of
#'   times it was selected divided by the total number of subsamples that were done 
#'   (total number of times the regression was performed).
#'
#'   We made use of the \code{stabs} package that implements lasso stability selection, 
#'   and adapted it to run randomized lasso stability selection. 
#'
#' @return  A \code{SummarizedExperiment} object where the rows are the 
#'   predictors and the columns are the different 
#'   regularization steps. The object consists of: \itemize{
#'      \item{assay}{:  \itemize{
#'        \item{selProb}{: the selection probabilities of the predictors 
#'        across the different regularization steps.}
#'        }
#'      }
#'      \item{metadata}{: list with: \itemize{
#'        \item{stabselParams}{: }
#'        \item{randStabselParams}{: }
#'        }
#'      }
#'   }
#'
#'@seealso \code{\link[stabs]{stabsel}}
#'
#'@references N. Meinshausen and P. Bühlmann (2010), Stability Selection, \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \strong{72}, 417–73. \cr
#'R.D. Shah and R.J. Samworth (2013), Variable Selection with Error Control: Another Look at Stability Selection, \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \strong{75}, 55–80. \cr
#'B. Hofner, L. Boccuto, and M. Göker (2015), Controlling False Discoveries in High-Dimensional Situations: Boosting with Stability Selection, \emph{BMC Bioinformatics}, \strong{16} 144.
#'
#' @importFrom stabs stabsel
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#'@export
randomized_stabsel <- function(x=NULL, y=NULL, weakness=0.8, cutoff=0.8, PFER=2, ...) {
  
    # checks
    if(!is(x, "matrix")) {
        stop("'x' must be a matrix")
    }
    .assertVector(y, type = "numeric")
    if(nrow(x) != length(y)) {
      stop("nrow of 'x' and length of 'y' are not equal. The rows of 
           x must be the same length and order as the elements in 'y'.")
    }
    if(!is.null(names(y)) && !is.null(rownames(x)) && names(y)!=rownames(x)) {
      stop("'x' and 'y' have different names. Make sure that the names are 
           identical and that the orders match.")
    }
    if(is.null(rownames(x))) {
        rownames(x) <- paste0("row_", 1:nrow(x))
    }
    if(is.null(colnames(x))) {
        colnames(x) <- paste0("col_", 1:ncol(x))
    }
  
    
    # run randomized lasso stability selection
    ss <- stabs::stabsel(x = x, y = y, fitfun = glmnet.randomized_lasso, 
                         args.fitfun = list(weakness = weakness), cutoff = cutoff, PFER = PFER, ...)

    
    # restructure as SummarizedExperiment object
    mdat <- list(stabselParams = list(cutoff = ss$cutoff, 
                                      selected = ss$selected, 
                                      max = ss$max, 
                                      q = ss$q, 
                                      PFER = ss$PFER, 
                                      specifiedPFER = ss$specifiedPFER, 
                                      p = ss$p, 
                                      B = ss$B, 
                                      sampling.type = ss$sampling.type, 
                                      assumption = ss$assumption, 
                                      call = ss$call), 
                 randStabselParams = list(weakness = weakness)
    )
    sel <- logical(length = ncol(x))
    sel[ss$selected] <- TRUE
    rdat <- DataFrame(sel)
    colnames(rdat) <- paste0("SelProbCutoff", cutoff) 
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(selProb = ss$phat), 
                                                     rowData = rdat, 
                                                     metadata = mdat)
   
    # return
    return(se)
    
}



