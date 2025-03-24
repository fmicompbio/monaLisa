#' @param allowNULL Logical, whether or not \code{NULL} is an acceptable
#'     value for \code{x}.
#'
#' @author Michael Stadler, Charlotte Soneson
#' @noRd
#' @keywords internal
#' @importFrom methods is
.assertScalar <- function(x,
                          type = NULL,
                          rngIncl = NULL,
                          rngExcl = NULL,
                          validValues = NULL,
                          allowNULL = FALSE) {

    .assertVector(x = x, type = type, rngIncl = rngIncl,
                  rngExcl = rngExcl, validValues = validValues,
                  len = 1, rngLen = NULL, allowNULL = allowNULL)

}

#' Utility function to check validity of vector variable values.
#'
#' This function provides a convenient way e.g. to check that provided
#' arguments to functions satisfy required criteria.
#'
#' @param x The variable to be checked
#' @param type The desired type of \code{x}.
#' @param rngIncl The allowed range of the (numeric) variable \code{x},
#'     including the endpoints.
#' @param rngExcl The allowed range of the (numeric) variable \code{x},
#'     excluding the endpoints.
#' @param validValues A vector with the allowed values of \code{x}.
#' @param len The required length of \code{x}.
#' @param rngLen The allowed range for the length of \code{x}.
#' @param allowNULL Logical, whether or not \code{NULL} is an acceptable
#'     value for \code{x}.
#'
#' @author Michael Stadler, Charlotte Soneson
#' @noRd
#' @keywords internal
#' @importFrom methods is
#' @importFrom cli cli_abort
.assertVector <- function(x,
                          type = NULL,
                          rngIncl = NULL,
                          rngExcl = NULL,
                          validValues = NULL,
                          len = NULL,
                          rngLen = NULL,
                          allowNULL = FALSE) {
    sc <- sys.calls()
    mycall <- sc[[length(sc)]]
    if (length(sc) >= 2 &&
        identical(as.character(sc[[length(sc) - 1]])[1], ".assertScalar")) {
        mycall <- sc[[length(sc) - 1]]
    }
    args <- lapply(mycall, \(x) if (is(x, "language")) deparse(x) else as.character(x))[-1]
    xname <- if ("x" %in% names(args)) args$x else "argument"

    ## Check arguments
    stopifnot(is.null(type) || (length(type) == 1L && is.character(type)))
    stopifnot(is.null(rngIncl) || (length(rngIncl) == 2L && is.numeric(rngIncl)))
    stopifnot(is.null(rngExcl) || (length(rngExcl) == 2L && is.numeric(rngExcl)))
    stopifnot(is.null(len) || (length(len) == 1L && is.numeric(len)))
    stopifnot(is.null(rngLen) || (length(rngLen) == 2L && is.numeric(rngLen)))
    stopifnot(is.logical(allowNULL) && length(allowNULL) == 1L)
    if (!is.null(rngIncl) && !is.null(rngExcl)) {
        cli_abort("{.arg rngIncl} and {.arg rngExcl} can not both be specified",
                  call = NULL)
    }

    ## If there are too many valid values, print only the first 15
    if (length(validValues) > 15) {
        vvPrint <- paste(c(validValues[seq_len(15)],
                           "...(truncated)"),
                         collapse = ", ")
    } else {
        vvPrint <- paste(validValues, collapse = ", ")
    }

    if (is.null(x)) {
        if (allowNULL) {
            return(invisible(TRUE))
        } else {
            cli_abort("{.arg {xname}} must not be {.code NULL}", call = NULL)
        }
    }

    if (is.null(type) && (!is.null(rngIncl) || !is.null(rngExcl))) {
        type <- "numeric"
    }

    if (!is.null(type) && !is(x, type)) {
        cli_abort("{.arg {xname}} must be of class {.cls {type}}", call = NULL)
    }

    if (!is.null(rngIncl)) {
        if (!is.null(validValues)) {
            if (any((x < rngIncl[1] | x > rngIncl[2]) & !(x %in% validValues))) {
                cli_abort(paste0(
                    "{.arg {xname}} must be between {rngIncl} (inclusive), ",
                    "or one of: {vvPrint}"), call = NULL)
            }
        } else {
            if (any(x < rngIncl[1] | x > rngIncl[2])) {
                cli_abort(paste0(
                    "{.arg {xname}} must be between {rngIncl} (inclusive)"),
                    call = NULL)
            }
        }
    } else if (!is.null(rngExcl)) {
        if (!is.null(validValues)) {
            if (any((x <= rngExcl[1] | x >= rngExcl[2]) & !(x %in% validValues))) {
                cli_abort(paste0(
                    "{.arg {xname}} must be between {rngExcl} (exclusive), ",
                    "or one of: {vvPrint}"), call = NULL)
            }
        } else {
            if (any(x <= rngExcl[1] | x >= rngExcl[2])) {
                cli_abort(paste0(
                    "{.arg {xname}} must be between {rngExcl} (exclusive)"),
                    call = NULL)
            }
        }
    } else {
        if (!is.null(validValues) && !all(x %in% validValues)) {
            cli_abort("All values in {.arg {xname}} must be one of: {vvPrint}",
                      call = NULL)
        }
    }


    if (!is.null(len) && length(x) != len) {
        cli_abort("{.arg {xname}} must have length {len}", call = NULL)
    }

    if (!is.null(rngLen) && (length(x) < rngLen[1] || length(x) > rngLen[2])) {
        cli_abort(paste0(
            "length of {.arg {xname}} must be between {rngLen} (inclusive)"),
            call = NULL)
    }

    return(invisible(TRUE))
}

#' Utility function that makes sure that packages are available
#' 
#' The function tries loading the namespaces of the packages given in
#' \code{pkgs}, and throws an exception with an informative error message if
#' that is not the case.
#' 
#' @param pkgs Character vector with package names. Can be either just a
#'   package name or a string of the form \code{"githubuser/packagename"} for
#'   packages hosted on GitHub.
#' @param suggestInstallation Logical scalar. If \code{TRUE}, include an
#'   expression to install the missing package(s) as part of the generated
#'   error message.
#' 
#' @author Michael Stadler, Charlotte Soneson
#' 
#' @noRd
#' @keywords internal
.assertPackagesAvailable <- function(pkgs, suggestInstallation = TRUE) {
  stopifnot(exprs = {
    is.character(pkgs)
    is.logical(suggestInstallation)
    length(suggestInstallation) == 1L
  })
  
  avail <- unlist(lapply(sub("^[^/]+/", "", pkgs),
                         function(pkg) {
                           requireNamespace(pkg, quietly = TRUE)
                         }))
  
  if (any(!avail)) {
    caller <- deparse(sys.calls()[[sys.nframe() - 1]])
    callerfunc <- sub("\\(.+$", "", caller)
    haveBioc <- requireNamespace("BiocManager", quietly = TRUE)
    msg <- paste0("The package", ifelse(sum(!avail) > 1, "s '", " '"),
                  paste(sub("^[^/]+/", "", pkgs[!avail]), collapse = "', '"),
                  "' ",
                  ifelse(sum(!avail) > 1, "are", "is"), " required for ",
                  callerfunc, "(), but not installed.\n")
    if (suggestInstallation) {
      msg <- paste0(msg,
                    "Install ", ifelse(sum(!avail) > 1, "them", "it"), " using:\n",
                    ifelse(haveBioc, "", "install.packages(\"BiocManager\")\n"),
                    "BiocManager::install(c(\"",
                    paste(pkgs[!avail], collapse = "\", \""), "\"))")
    }
    stop(msg, call. = FALSE)
  }
  
  invisible(TRUE)
}

#' Generate console output messages
#'
#' This is a drop-in replacement for \code{base::message}, which only
#' creates a message if \code{verbose} exists in the calling environment and
#' is set to \code{TRUE}. It also supports inline-markup via the \code{cli}
#' package.
#'
#' @param message The message to be written to the console. It will be
#'     forwarded to \code{\link[cli]{cli_progress_step}} and thus supports
#'     inline markup (see \code{\link[cli]{inline-markup}}).
#' @param noTimer Logical scalar. If \code{FALSE} (the default), the message is
#'     generated using \code{\link[cli]{cli_progress_step}}, which will first
#'     show it with an "info" icon and then again with a "check" icon and timing
#'     information when the next message is generated or the function terminates.
#'     If \code{TRUE}, the message is generated using
#'     \code{\link[cli]{cli_alert_info}}, which means it will be show only once
#'     with an "info" icon.
#' @param ... Additional arguments passed to \code{\link[cli]{cli_progress_step}}
#'     (ignored if \code{noTimer = TRUE}).
#'
#' @importFrom cli cli_progress_step cli_alert_info
#'
#' @noRd
#' @keywords internal
.message <- function(message, noTimer = FALSE, ...) {
    # Try to get 'verbose' from the calling environment
    env <- parent.frame()
    verbose <- tryCatch(get("verbose", envir = env),
                        error = function(e) FALSE)
    if (verbose) {
        if (noTimer) {
            cli_alert_info(text = message, .envir = env)
        } else {
            cli_progress_step(msg = message, .envir = env, ...)
        }
    }
}

#' Utility function to check validity of R colors.
#'
#' This function checks if the provided characters (e.g. arguments to other
#' functions) are valid colors in R.
#'    
#' @param x character vector of colors.    
#' @param ... additional arguments to \code{.assertVector}
#'    
#' @importFrom grDevices col2rgb
#'
#' @noRd
#' @keywords internal
.assertColor <- function(x, ...) {
    .assertVector(x = x, type = "character", ...)  
    
    isColor <- function(char) {
        out <- tryCatch(col2rgb(char), error = function(e) {
            FALSE
        })
        if (is.matrix(out)){
            out <- TRUE
        }
        out
    }  
    
    out <- vapply(X = x, FUN = isColor, FUN.VALUE = logical(1))
    if (any(!out)) {
        cli_abort("The following is not a valid color in R: {.emph {x[!out]}}. Please
              provide a string that corresponds to a valid color.")
    } 
    
    return(invisible(TRUE))
}

