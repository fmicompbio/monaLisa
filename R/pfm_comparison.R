# compare two PFMs of any length (padd left/right with background positions)
# score := correlation of single base frequencies of aligned and padded matrices
# (internal function used by motifSimilarity)
#' @importFrom stats cor
#' 
#' @keywords internal
.compareMotifPair <- function(m1, m2) {
    bestScore <- -2
    bestOffset <- 0
    bestDirection <- ""

    ## reverse complement
    rv2 <- rev(m2)
    attributes(rv2) <- attributes(m2)

    len1 <- ncol(m1)
    len2 <- ncol(m2)

    ## offset := start(m1) - start(m2)
    ##        == left-padding of m1 (if > 0)
    ##        == -1 * left-padding of m2 (if < 0)
    for (offset in seq(-len1 + 1L, len2 - 1L)) { # minimal overlap of 1
        # padding of matrices
        mm1 <- cbind(matrix(0.25, nrow = 4, ncol = max(0, offset)),
                     m1,
                     matrix(0.25, nrow = 4, 
                            ncol = max(0, len2 - (len1 + offset))))
        mm2 <- cbind(matrix(0.25, nrow = 4, ncol = max(0, -offset)),
                     m2,
                     matrix(0.25, nrow = 4, 
                            ncol = max(0, (len1 + offset) - len2)))
        mm2r <- cbind(matrix(0.25, nrow = 4, ncol = max(0, -offset)),
                      rv2,
                      matrix(0.25, nrow = 4, 
                             ncol = max(0, (len1 + offset) - len2)))

        score <- cor(as.vector(mm1), as.vector(mm2))
        rvScore <- cor(as.vector(mm1), as.vector(mm2r))

        if (rvScore > bestScore) {
            bestScore <- rvScore
            bestOffset <- offset
            bestDirection <- "revcomp"
        }
        if (score > bestScore) {
            bestScore <- score
            bestOffset <- offset
            bestDirection <- "forward"
        }
    }

    return(list(bestScore = bestScore, bestOffset = bestOffset, 
                bestDirection = bestDirection))
}

# compare a PFM to all k-mer of any length (padd left/right with 
# background positions)
# score := maximal probability of observing k-mer under (potentially padded) PFM
# (internal function used by motifKmerSimilarity)
#' 
#' @keywords internal
.compareMotifKmer <- function(m, kmers) {
    bestScore <- rep(-2, length(kmers))
    bestOffset <- rep(0, length(kmers))

    len1 <- ncol(m)
    len2 <- nchar(kmers[1])
    j <- seq.int(len2)

    ## transform k-mers into numeric indices
    kmersL <- strsplit(x = kmers, split = "", fixed = TRUE)
    kmersN <- utils::relist(match(unlist(kmersL), c("A","C","G","T")), kmersL)
    ## offset := start(m) - start(kmers[i])
    ##        == left-padding of m (if > 0)
    ##        == -1 * left-padding of kmers[i] (if < 0)
    for (offset in seq(-len1 + 1L, len2 - 1L)) {
        # minimal overlap of 1
        # padding of matrix and kmers
        mm <- cbind(matrix(0.25, nrow = 4, ncol = max(0, offset)),
                    m[, seq.int(min(len1 + offset, 
                                    len2 - offset, len2)) + max(0, -offset)],
                    matrix(0.25, nrow = 4, 
                           ncol = max(0, len2 - (len1 + offset))))
        score <- unlist(lapply(kmersN, function(i) prod(mm[cbind(i, j)])))
        b <- score > bestScore
        if (any(b)) {
            bestScore[b] <- score[b]
            bestOffset[b] <- offset
        }
    }

    return(list(bestScore = bestScore, bestOffset = bestOffset))
}

#' @title Calculate similarities between pairs of motifs.
#'
#' @description For each pair of motifs, calculate the similarity defined as the
#'   maximal Pearson's correlation coefficient between base frequencies over all
#'   possible shifts (relative positions of the two matrices with at least one
#'   overlapping position). If necessary matrices are padded on the sides with
#'   background base frequencies (assuming all bases to have a frequency of
#'   0.25) to enable comparison of all positions in both matrices.
#'
#' @param x Either a \code{\link[TFBSTools]{PFMatrixList}}, or a character
#'   scalar with a file containing motifs in HOMER format (used directly
#'   \code{method = "HOMER"}, loaded into a
#'   \code{\link[TFBSTools]{PFMatrixList}} by \code{\link{homerToPFMatrixList}}
#'   for \code{method = "R"}).
#' @param y Either a \code{\link[TFBSTools]{PFMatrixList}} or \code{NULL}
#'   (default). If \code{y = NULL}, then similarities will be calucalted for all
#'   pairs of motifs within \code{x}. Otherwise, \code{method} must be
#'   \code{"R"} and similarities will be calculated between any motif from
#'   \code{x} to any motif from \code{y}.
#' @param method A character scalar specifying the method for similarity
#'   calculations. Either \code{"R"} (pure R implementation) or \code{"HOMER"}
#'   (will call the \code{compareMotifs.pl} script from HOMER). Results are
#'   identical (apart from rounding errors), and the R implementation is
#'   usually faster and can be parallelized (\code{BPPARAM} argument).
#' @param homerfile Path to the HOMER script \code{compareMotifs.pl} (only used
#'   for \code{method = "HOMER"}.
#' @param homerOutfile A character scalar giving the file to save the similarity
#'   scores (only for \code{metho = "HOMER"}). If \code{NULL}, scores will be
#'   stored into a temporary file.
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'     instance determining the parallel back-end to be used during evaluation
#'     (only used for \code{method = "R"}).
#' @param verbose A logical scalar. If \code{TRUE}, report on progress.
#'
#' @return A matrix of Pearson's correlation coefficients for each pair of
#'   motifs.
#'
#' @examples 
#' m <- rbind(A = c(12,  0,  0),
#'            C = c( 3,  2,  0),
#'            G = c( 0, 14,  0),
#'            T = c( 0,  0, 15))
#' pfms <- TFBSTools::PFMatrixList(
#'     TFBSTools::PFMatrix(name = "m1", profileMatrix = m),
#'     TFBSTools::PFMatrix(name = "m2", profileMatrix = m + 10),
#'     TFBSTools::PFMatrix(name = "m3", profileMatrix = m[, 3:1])
#' )
#' motifSimilarity(pfms)
#'
#' @seealso \code{\link[BiocParallel]{bplapply}} used for parallelization for
#'   \code{method = "R"},
#'   documentation of HOMER's \code{compareMotifs.pl} for details on
#'   \code{method = "HOMER"}.
#'
#' @importFrom BiocParallel bplapply SerialParam bpnworkers
#' 
#' @export
motifSimilarity <- function(x, y = NULL, method = c("R", "HOMER"),
                            homerfile = findHomer("compareMotifs.pl"),
                            homerOutfile = NULL, BPPARAM = SerialParam(),
                            verbose = FALSE) {
    ## branch by method
    stopifnot(exprs = {is.logical(verbose); length(verbose) == 1L})
    method <- match.arg(method)
    if (method == "R") {
        ## pre-flight checks for "R"
        if (is.character(x) && length(x) == 1L && file.exists(x)) {
            if (verbose) {
                message("reading motifs from ", basename(x))
            }
            x <- homerToPFMatrixList(x)
        }
        stopifnot(exprs = {
            is(x, "PFMatrixList")
            is.null(y) || is(y, "PFMatrixList")
            is(BPPARAM, "BiocParallelParam")
        })

        if (bpnworkers(BPPARAM) > 2 && is.null(y)) {
            y <- x
        }

        xm <- lapply(TFBSTools::Matrix(x), 
                     function(x) sweep(x, 2, colSums(x), "/"))

        if (is.null(y)) { # compare x to itself, n*(n-1)/2 comparisons
            if (verbose) {
                message("calculating ", length(xm) * (length(xm) - 1) / 2,
                        " similarities...", appendLF = FALSE)
            }
            M <- matrix(NA, nrow = length(xm), ncol = length(xm), 
                        dimnames = list(name(x), name(x)))
            diag(M) <- 1.0
            for (i in seq.int(length(x) - 1L)) {
                for (j in seq(i + 1, length(x))) {
                    M[i, j] <- M[j, i] <- 
                        .compareMotifPair(xm[[i]], xm[[j]])$bestScore
                }
            }
            if (verbose) {
                message("done")
            }
        } else {#           compare x to y, n*m comparisons
            ym <- lapply(TFBSTools::Matrix(y), 
                         function(x) sweep(x, 2, colSums(x), "/"))
            if (verbose) {
                message(
                    "calculating ", length(xm) * length(ym),
                    " similarities using ", bpnworkers(BPPARAM),
                    if (bpnworkers(BPPARAM) > 1) " cores..." else " core...",
                    appendLF = FALSE
                )
            }
            if (bpnworkers(BPPARAM) > 1) {
                M <- do.call(rbind, bplapply(seq_along(xm), function(i) {
                    unlist(lapply(
                        seq_along(ym), 
                        function(j) .compareMotifPair(xm[[i]], 
                                                      ym[[j]])$bestScore
                    ))
                }, BPPARAM = BPPARAM))
                dimnames(M) <- list(name(x), name(y))
            } else {
                M <- matrix(NA, nrow = length(xm), ncol = length(ym), 
                            dimnames = list(name(x), name(y)))
                for (i in seq_along(xm)) {
                    for (j in seq_along(ym)) {
                        M[i, j] <- .compareMotifPair(xm[[i]], ym[[j]])$bestScore
                    }
                }
            }
            if (verbose) {
                message("done")
            }
        }

    } else if (method == "HOMER") {
        ## pre-flight checks for "HOMER"
        stopifnot(exprs = {is.character(x); length(x) == 1L; file.exists(x)})
        stopifnot(exprs = {!is.na(homerfile); is.character(homerfile); 
            length(homerfile) == 1L; file.exists(homerfile)})
        if (is.null(homerOutfile)) {
            homerOutfile <- tempfile(fileext = ".simmat")
        }
        stopifnot(exprs = {is.character(homerOutfile); 
            length(homerOutfile) == 1L; !file.exists(homerOutfile)})

        ## run
        if (verbose) {
            message("running compareMotifs.pl...")
        }
        system2(command = homerfile,
                args = sprintf("%s test -matrix %s", x, homerOutfile),
                stdout = FALSE, stderr = FALSE)
        M <- as.matrix(read.delim(homerOutfile, row.names = 1))
    }
    return(M)
}


#' @title Calculate similarities between motifs and k-mers.
#'
#' @description For each motif, calculate it's similarity to all k-mers of
#'   length \code{kmerLen}, defined as the maximal probability of observing the
#'   k-mer given the base frequencies of the motif (the maximum is taken over
#'   for all possible ungapped alignments between motif and k-mer). If necessary
#'   matrices are padded on the sides with background base frequencies (assuming
#'   all bases to have a frequency of 0.25).
#'
#' @param x Either a \code{\link[TFBSTools]{PFMatrixList}}, or a character
#'   scalar with a file containing motifs in HOMER format (used directly
#'   \code{method = "HOMER"}, loaded into a
#'   \code{\link[TFBSTools]{PFMatrixList}} by \code{\link{homerToPFMatrixList}}
#'   for \code{method = "R"}).
#' @param kmerLen A \code{numeric} scalar giving the k-mer length.
#' @param kmers Either a character vector of k-mers for which to calculate 
#'   the similarity to each motif, or \code{NULL}, in which case all k-mers
#'   of length \code{kmerLen} are used.
#' @param includeRevComp A \code{logical} scalar. If set to \code{TRUE}, each 
#'   k-mer as well as its reverse complement is compared to each motif, and the 
#'   larger of the two similarities is returned. 
#' @param BPPARAM An optional \code{\link[BiocParallel]{BiocParallelParam}}
#'   instance determining the parallel back-end to be used during evaluation.
#' @param verbose A logical scalar. If \code{TRUE}, report on progress.
#'
#' @return A matrix of probabilties for each motif - k-mer pair.
#' 
#' @examples 
#' m <- rbind(A = c(12,  0,  0),
#'            C = c( 3,  2,  0),
#'            G = c( 0, 14,  0),
#'            T = c( 0,  0, 15))
#' pfms <- TFBSTools::PFMatrixList(
#'     TFBSTools::PFMatrix(name = "m1", profileMatrix = m),
#'     TFBSTools::PFMatrix(name = "m2", profileMatrix = m[, 3:1])
#' )
#' motifKmerSimilarity(pfms, kmerLen = 3)[, c("AGT", "TGA")]
#'
#' @seealso \code{\link[BiocParallel]{bplapply}} used for parallelization.
#'
#' @importFrom BiocParallel bplapply SerialParam bpnworkers
#' @importFrom Biostrings DNA_ALPHABET
#' 
#' @export
motifKmerSimilarity <- function(x,
                                kmerLen = 5,
                                kmers = NULL,
                                includeRevComp = FALSE, 
                                BPPARAM = SerialParam(),
                                verbose = FALSE) {
    ## pre-flight checks
    .assertScalar(x = verbose, type = "logical")
    .assertScalar(x = includeRevComp, type = "logical")
    if (is.character(x) && length(x) == 1L && file.exists(x)) {
        if (verbose) {
            message("reading motifs from ", basename(x))
        }
        x <- homerToPFMatrixList(x)
    }
    .assertScalar(x = kmerLen, type = "numeric", rngExcl = c(0, Inf))
    stopifnot(exprs = {
        is.null(kmers) || 
            (.assertVector(kmers, type = "character") && 
                 all(unlist(strsplit(kmers, "")) %in% Biostrings::DNA_ALPHABET))
    })
    stopifnot(exprs = {
        is(x, "PFMatrixList")
        round(kmerLen, 0L) == kmerLen
        is(BPPARAM, "BiocParallelParam")
    })

    xm <- lapply(TFBSTools::Matrix(x), function(x) sweep(x, 2, colSums(x), "/"))
    if (is.null(kmers)) {
        kmers <- Biostrings::mkAllStrings(c("A", "C", "G", "T"), kmerLen)
    }

    if (verbose) {
        message("calculating ", length(xm) * length(kmers),
                " similarities using ", bpnworkers(BPPARAM),
                if (bpnworkers(BPPARAM) > 1) " cores..." else " core...",
                appendLF = FALSE)
    }
    
    if (includeRevComp) {
        kmersrevcomp <- as.character(Biostrings::reverseComplement(
            Biostrings::DNAStringSet(kmers)
        ))
        kmersall <- union(kmers, kmersrevcomp)
    } else {
        kmersall <- kmers
    }
    
    M0 <- do.call(rbind, bplapply(
        xm, function(m) .compareMotifKmer(m = m, kmers = kmersall)$bestScore,
        BPPARAM = BPPARAM))
    dimnames(M0) <- list(name(x), kmersall)
    
    ## Get the matrix for the originally specified k-mers
    M <- M0[, kmers, drop = FALSE]
    
    if (includeRevComp) {
        ## Get the matrix for the corresponding reverse complement k-mers, 
        ## and take the maximum similarity for each motif
        Mrev <- M0[, kmersrevcomp, drop = FALSE]
        M <- pmax(M, Mrev)
    }
    
    if (verbose) {
        message("done")
    }

    return(M)
}

