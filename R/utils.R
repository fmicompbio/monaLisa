# utility functions used to check parameter conformity
# - 'x' is the argument to be checked
# - additional arguments define conditions
# - a violated condition will throw an exception using stop()
# - if no condition is violated, TRUE is returned (invisibly)

.assertScalar <- function(x,
                          type = NULL,
                          rngIncl = NULL,
                          rngExcl = NULL) {
    args <- lapply(sys.call()[-1], as.character)
    xname <- if ("x" %in% names(args)) args$x else "argument"

    if (length(x) != 1L) {
        stop("'", xname, "' must be a scalar value (length one)")
    }
    
    if (is.null(type) && (!is.null(rngIncl) || !is.null(rngExcl))) {
        type <- "numeric"
    }

    if (!is.null(type) && !is(x, type)) {
        stop("'", xname, "' must be of type '", type, "'")
    }
    
    if (!is.null(rngIncl) && is.numeric(rngIncl) && length(rngIncl) == 2L &&
        (x < rngIncl[1] || x > rngIncl[2])) {
        stop("'", xname, "' must be within [", rngIncl[1], ",", rngIncl[2], "] (inclusive)")
    }
    
    if (!is.null(rngExcl) && is.numeric(rngExcl) && length(rngExcl) == 2L &&
        (x <= rngExcl[1] || x >= rngExcl[2])) {
        stop("'", xname, "' must be within (", rngExcl[1], ",", rngExcl[2], ") (exclusive)")
    }
    
    return(invisible(TRUE))
}


.assertVector <- function(x,
                          type = NULL,
                          rngIncl = NULL,
                          rngExcl = NULL,
                          len = NULL) {
    args <- lapply(sys.call()[-1], as.character)
    xname <- if ("x" %in% names(args)) args$x else "argument"
    
    if (is.null(type) && (!is.null(rngIncl) || !is.null(rngExcl))) {
        type <- "numeric"
    }
    
    if (!is.null(type) && !is(x, type)) {
        stop("'", xname, "' must be of type '", type, "'")
    }
    
    if (!is.null(rngIncl) && is.numeric(rngIncl) && length(rngIncl) == 2L &&
        any(x < rngIncl[1] | x > rngIncl[2])) {
        stop("values in '", xname, "' must be within [", rngIncl[1], ",", rngIncl[2], "] (inclusive)")
    }
    
    if (!is.null(rngExcl) && is.numeric(rngExcl) && length(rngExcl) == 2L &&
        any(x <= rngExcl[1] | x >= rngExcl[2])) {
        stop("values in '", xname, "' must be within (", rngExcl[1], ",", rngExcl[2], ") (exclusive)")
    }

    if (!is.null(len) && is.numeric(len) && length(len) == 1L && length(x) != len) {
        stop("'", xname, "' must have length ", len)
    }
    
    return(invisible(TRUE))
}



