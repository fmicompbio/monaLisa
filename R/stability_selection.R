#' @importFrom stabs stabsel
#' @importFrom glmnet glmnet
#' @importFrom TFBSTools PWMatrixList
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
#'@author Dania Machlab
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
#'@description This function runs randomized lasso stability selection as presented by Meinshausen
#'and Bühlmann (2010) and with the improved error bounds introduced by Shah and Samworth (2013). The function
#'uses the \code{stabsel} function from the \code{stabs} package, but implements the randomized lasso version.
#'
#'@param x the predictor matrix.
#'@param y the response vector.
#'@param weakness parameter that can take a value between 0 and 1. It affects how strict the method will be in
#'selecting predictors. It is set to 0.8 by default. The closer it goes to 0 the more stringent the selection. A
#'weakness value of 1 is identical to performing lasso stability selection (not the randomized version).
#'@param cutoff value between 0 and 1 which is the cutoff for the selection probability. Any variable
#'with a selection probability that is higher than the set cutoff will be selected. It is set to 0.8 by default.
#'@param PFER the absolute number of false positives that we allow for in the final list of
#'selected variables. For details see Meinshausen and Bühlmann (2010).
#'@param ... additional parameters for \code{glmnet}.
#'
#'@return a \code{stabsel} object. It contains the following:
#'\itemize{
#'  \item phat - a matrix containing selection probabilities, where columns are regularization steps and rows
#'  are predictors.
#'  \item selected - a vector containing the selected predictors.
#'  \item max - maximum of selection probabilities.
#'  \item cutoff - the selection probability cutoff.
#'  \item q - average number of selected variables used.
#'  \item PFER - the realized upper bound for the per-family error rate (number of falsely selected predictors among the group of selected predictors).
#'  \item specifiedPFER - the set upper bound for the per-family error rate (number of falsely selected predictors among the group of selected predictors).
#'  \item p - the total number of predictors.
#'  \item B - the number of subsamples.
#'  \item sampling.type - sampling type used for stability selection.
#'  \item assumption - assumptions made on the selection probabilities.
#'  \item call - the call.
#'}
#'
#'@details Randomized lasso stability selection runs randomized lasso several times on subsamples of
#'the response variable and predictor matrix (stability selection). N/2 elements from the response variable are
#'randomly chosen, where N is the length of the vector. Their corresponsing section of the predictor matrix is
#'also chosen, and the \code{glmnet.randomized_lasso} function is applied. This is done multiple times, and
#'results in selection probabilities for each predictor. The probability of a specific predictor is the number of
#'times it was selected divided by the total number of subsamples that were done (total number of times the
#'regression was performed).
#'
#'We make use of the \code{stabs} package that implements lasso stability selection, and adapt it to run the
#'randomized lasso stability selection. The output is an object of type \code{stabsel}.
#'
#'@seealso \code{\link[stabs]{stabsel}}
#'
#'@references N. Meinshausen and P. Bühlmann (2010), Stability Selection, \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \strong{72}, 417–73. \cr
#'R.D. Shah and R.J. Samworth (2013), Variable Selection with Error Control: Another Look at Stability Selection, \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \strong{75}, 55–80. \cr
#'B. Hofner, L. Boccuto, and M. Göker (2015), Controlling False Discoveries in High-Dimensional Situations: Boosting with Stability Selection, \emph{BMC Bioinformatics}, \strong{16} 144.
#'
#'
#'@author Dania Machlab
#'@export
randomized_stabsel <- function(x=x, y=y, weakness=0.8, cutoff=0.8, PFER=2, ...) {
  stabs::stabsel(x = x, y = y, fitfun = glmnet.randomized_lasso, args.fitfun = list(weakness = weakness), cutoff = cutoff, PFER = PFER, ...)
}




#' #'@title TFBS Matrix
#' #'
#' #'@description This function takes in a \code{GRanges} that is the output of \code{findMotifHits} and returns a
#' #'matrix contaiting the number of transcription factor bisnding sites (TFBS) per motif across a given set of genomic ranges.
#' #'
#' #'@param TFBS_gr the output of \code{findMotifHits} contaiting the TF binding locations across specified genomic regions.
#' #'@param subject_gr a \code{GRanges} object showing the positions of the geominc regions that have been scanned for the TFs.
#' #'This corresponds to the \code{GRanges} \code{subject} parameter used in the \code{findMotifHits} function.
#' #'@param PWMs \code{PWMatrixList} or \code{PWMatrix} object of the used TFs. This corresponds to the \code{query}
#' #'parameter used in the \code{findMotifHits} function.
#' #'@param Ncpu number of CPUs to use (set to 1 by default).
#' #'
#' #'@return a matrix containing the number of binding sites each TF has across the genomic regions.
#' #'
#' #'@author Dania Machlab
#' #'@export
#' get_numberOfTFBS_perSeqName <- function(TFBS_gr, subject_gr, PWMs, Ncpu=1L) {
#' 
#'   # TODO make sure PWM names and peak names are  unique
#'   # motif names instead of gene?
#' 
#'   ## checks
#'   stopifnot(base::inherits(TFBS_gr, "GRanges"))
#'   stopifnot(all(colnames(as.data.frame(TFBS_gr)) == c("seqnames", "start", "end", "width", "strand", "matchedSeq", "pwmname", "score")))
#'   stopifnot(base::inherits(PWMs, what="PWMatrixList") | base::inherits(PWMs, what="PWMatrix"))
#'   if (base::inherits(PWMs, what="PWMatrix")) {PWMs <- TFBSTools::PWMatrixList(PWMs)}
#'   # if(!all(TFBS_gr$pwmname %in% sapply(PWMs, function(x){name(x)}))) {stop("PWMs missing motifs found in TFBS_gr")}
#' 
#'   ## for each motif count number of TFBS per seqName
#'   seqs <- as.character(seqnames(TFBS_gr))
#'   TFs <- as.character(TFBS_gr$pwmname)
#'   s <- split(TFs, seqs)
#'   l <- lapply(s, function(x){table(x)})
#' 
#'   ## output full matrix
#'   if (is.null(names(subject_gr))) {
#'     names(subject_gr) <- paste0("row_", seq(from = 1, to = length(subject_gr), by = 1))
#'   }
#' 
#'   ## rbind vectors
#'   m <- do.call(rbind, parallel::mclapply(mc.cores = Ncpu, X = l, FUN = function(x) {
#'     full_motif_vec <- numeric(length(PWMs))
#'     names(full_motif_vec) <- sapply(PWMs, function(x){name(x)})
#'     df <- as.data.frame(x)
#'     motifs <- as.character(df$x)
#'     full_motif_vec[motifs] <- df$Freq
#'     full_motif_vec
#'   }))
#' 
#'   ## order to match subject_gr
#'   subject_peaks <- names(subject_gr)[names(subject_gr) %in% rownames(m)] # remove peaks that have 0 TFBS in all columns
#'   o <- match(subject_peaks, rownames(m))
#'   m <- m[o, ]
#' 
#'   ## return matrix
#'   m
#' 
#' }







