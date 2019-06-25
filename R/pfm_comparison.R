# comparison of PFMs, based on Homer's compareMotifs.pl

# compare two PFMs of any length (padd left/right with background positions)
compareMotifPair <- function(m1, m2) {
    stopifnot(is.matrix(m1) && is.matrix(m2) &&
                  nrow(m1) == 4L && nrow(m2) == 4L &&
                  all.equal(rep(1.0, ncol(m1)), colSums(m1)) &&
                  all.equal(rep(1.0, ncol(m2)), colSums(m2)))

    bestScore <- -1e10
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
                     matrix(0.25, nrow = 4, ncol = max(0, len2 - (len1 + offset))))
        mm2 <- cbind(matrix(0.25, nrow = 4, ncol = max(0, -offset)),
                     m2,
                     matrix(0.25, nrow = 4, ncol = max(0, (len1 + offset) - len2)))
        mm2r <- cbind(matrix(0.25, nrow = 4, ncol = max(0, -offset)),
                      rv2,
                      matrix(0.25, nrow = 4, ncol = max(0, (len1 + offset) - len2)))

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

    return(list(bestScore = bestScore, bestOffset = bestOffset, bestDirection = bestDirection))
}

# compare all pairs of motifs between two sets of PFMs
compareAllMotifPairs <- function(x, y = NULL, Ncpu = 1L) {
    stopifnot(is(x, "PFMatrixList") && (is.null(y) || is(y, "PFMatrixList")))
    stopifnot(is.numeric(Ncpu) && length(Ncpu) == 1 && Ncpu > 0)

    xm <- lapply(TFBSTools::Matrix(x), function(x) sweep(x, 2, colSums(x), "/"))

    if (Ncpu > 2 && is.null(y)) {
        y <- x
    }

    if (is.null(y)) { # compare x to itself, n*(n-1)/2 comparisons
        M <- matrix(NA, nrow = length(xm), ncol = length(xm), dimnames = list(name(x), name(x)))
        diag(M) <- 1.0
        for (i in seq.int(length(x) - 1L)) {
            for (j in seq(i + 1, length(x))) {
                M[i, j] <- M[j, i] <- compareMotifPair(xm[[i]], xm[[j]])$bestScore
            }
        }
    } else {         # compare x to y, n*m comparisons
        ym <- lapply(TFBSTools::Matrix(y), function(x) sweep(x, 2, colSums(x), "/"))
        if (Ncpu > 1) {
            M <- do.call(rbind, parallel::mclapply(seq_along(xm), function(i) {
                unlist(lapply(seq_along(ym), function(j) compareMotifPair(xm[[i]], ym[[j]])$bestScore))
            }, mc.cores = Ncpu))
            dimnames(M) <- list(name(x), name(y))
        } else {
            M <- matrix(NA, nrow = length(xm), ncol = length(ym), dimnames = list(name(x), name(y)))
            for (i in seq_along(xm)) {
                for (j in seq_along(ym)) {
                    M[i, j] <- compareMotifPair(xm[[i]], ym[[j]])$bestScore
                }
            }
        }
    }

    return(M)
}



#' @title Calculate similarity matrix of motifs.
#'
#' @description Run the HOMER script compareMotifs.pl (with default options) to get a similarity matrix
#'     of all motifs. For details, see the HOMER documentation.
#'
#'
#' @param motifFile A file with HOMER formatted PWMs as input for compareMotifs.pl.
#' @param homerdir Path to the HOMER binary directory.
#' @param outfile A file to save the similarity scores.
#' @return A matrix of Pearson correlations for each pairwise comparison of
#'     motifs in motifFile.
#'
#' @export
motifSimilarity <- function(motifFile, homerdir, outfile){

    stopifnot(is.character(motifFile) && length(motifFile) == 1L && file.exists(motifFile))
    stopifnot(is.character(homerdir) && length(homerdir) == 1L && file.exists(homerdir))
    homerfile = findHomer("compareMotifs.pl", dirs = homerdir)
    #run
    message("running compareMotifs.pl...")
    system(sprintf("%s %s test -matrix %s", homerfile, motifFile, outfile), intern=TRUE)
    #/work/gbioinfo/Appz/Homer/Homer-4.8/bin/compareMotifs.pl test.motif test -matrix testmat
    as.matrix(read.delim(outfile, row.names=1))
}

# pfms <- TFBSTools::getMatrixSet(JASPAR2018, list(tax_group = "vertebrates"))
# name(pfms[c("MA0139.1", "MA1102.1", "MA1104.1")])
# tmp <- lapply(Matrix(pfms[c("MA0139.1", "MA1102.1", "MA1104.1")]), function(x) sweep(x, 2, colSums(x), "/"))
# for (i in 1:2) { for (j in i:3) { print(paste(i,"-",j,":", lisa:::compareMotifPair(tmp[[i]], tmp[[j]])$bestScore)) } }
# compareAllMotifPairs(pfms[c("MA0139.1", "MA1102.1", "MA1104.1")])
#
# tmp2 <- lapply(list(Matrix(pfms[[1]]), cbind(Matrix(pfms[[1]])[4:1,][,6:1],matrix(1,nrow=4, ncol=6))), function(x) sweep(x, 2, colSums(x), "/"))
# compareMotifPair(tmp2[[1]], tmp2[[2]])
#
# length(pfms) # 579
#
# motiffile <- tempfile(fileext = ".motif")
# dumpJaspar(motiffile, pkg = "JASPAR2018")
# system.time(M0 <- clusterPWMs(motifFile = motiffile,
#             homerdir ="/work/gbioinfo/Appz/Homer/Homer-4.10.4/bin/",
#             outfile = tempfile(fileext = ".simmat")))
#    user  system elapsed
# 568.937   0.313 570.500
#
# system.time(M <- compareAllMotifPairs(pfms))
#    user  system elapsed
# 228.340   4.810 233.651
#
# all.equal.numeric(M0, M, check.attributes = FALSE) # "Mean relative difference: 0.0006257828"
# summary(diag(cor(M0, M)))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.9991  1.0000  1.0000  1.0000  1.0000  1.0000
#
# system.time(M2 <- compareAllMotifPairs(pfms, Ncpu = 30))
#    user  system elapsed
# 471.356  71.179  23.236
# all.equal(M, M2) # TRUE
#
# library(microbenchmark)
# microbenchmark(lisa:::compareMotifPair(tmp[[1]], tmp[[2]]))

# SimMat[c(51,52,166),c(51,52,166)]
# CTCF_P49711_ChIPseq CTCFL_Q8NI51_ChIPseq GATA6_Q92908_ChIPseq
# CTCF_P49711_ChIPseq            1.0000000            0.8431996            0.1991379
# CTCFL_Q8NI51_ChIPseq           0.8431996            1.0000000            0.1426231
# GATA6_Q92908_ChIPseq           0.1991379            0.1426231            1.0000000
