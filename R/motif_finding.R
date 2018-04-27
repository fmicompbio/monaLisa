#' @importFrom Biostrings matchPWM DNA_BASES DNAStringSet
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom TFBSTools Matrix PWMatrix PWMatrixList name
#' @importFrom methods .valueClassTest
NULL

# internal function:  dump PWMatrixList to file in homer2 format
.dumpPWMsToHomer2File <- function(pwmL, fname, relscore=0.8, absscore=NULL) {
    stopifnot(inherits(pwmL, "PWMatrixList"))
    stopifnot(is.null(absscore) || is.numeric(absscore))
    if (is.numeric(absscore) && length(absscore) == 1L)
        absscore <- rep(absscore, length(pwmL))
    fh <- file(fname, "wb")
    for (i in seq_along(pwmL)) {
        pwm <- TFBSTools::Matrix(pwmL[[i]])
        scorecut <- if (is.null(absscore)) relscore * sum(log(apply(pwm, 2, max) / 0.25)) else absscore[i]
        cat(sprintf(">%s\t%s\t%.2f\n",
                    paste(apply(pwm, 2, function(x) { rownames(pwm)[which.max(x)] }), collapse = ""),
                    TFBSTools::name(pwmL[[i]]), scorecut),  file = fh, append = TRUE)
        write.table(file = fh, t(pwm), row.names = FALSE, col.names = FALSE,
                    sep = "\t", quote = FALSE, append = TRUE)
    }
    close(fh)
    return(fname)
}

#' @title Find motif matches in sequences.
#'
#' @description \code{findMotifHits} scans sequences (either provided
#' as a file, an R object or genomic coordinates)  for matches to
#' positional weight matrices (provided as a file or as R objects)
#'
#' @param query The motifs to search for, either a
#'     \itemize{
#'         \item{\code{character(1)}}{ with the path and file name of a motif
#'               file with PWM in HOMER format (currently only
#'               supported for \code{method="homer2"})}
#'         \item{\code{PWMatrix}}{ with a single PWM}
#'         \item{\code{PWMatrixList}}{ with several PWMs to search for.}
#'      }
#'
#' @param subject The sequences to be searched, either a
#'     \itemize{
#'         \item{\code{character}}{ with the path and file name of a sequence
#'               file with DNA sequences in FASTA format}
#'         \item{\code{DNAString}}{ with a single sequence}
#'         \item{\code{DNAStringSet}}{ with several sequences}
#'         \item{\code{GRanges}}{ object with the genomic coordinates
#'             of the sequences to be searched.}
#'     }
#'
#' @param min.score The minimum score for counting a match. Can be given as
#'     a character string containing a percentage (e.g. "85%") of the highest
#'     possible score or as a single number.
#'
#' @param method The internal method to use for motif searching. One of
#'     \itemize{
#'         \item{\code{homer2}}{call to the homer2 binary}
#'         \item{\code{"matchPWM"}}{using Biostrings::matchPWM}
#'     }
#'
#' @param homerfile Path and file name of the \code{homer2} binary.
#'
#' @param Ncpu Number of parallel threads (only for \code{method="homer2"}, default: 1).
#'
#' @param ... Additional arguments for specific methods.
#'
#' @return A \code{GRanges} object with the matches to \code{query} in \code{subject}.
#'
#' @seealso \code{\link[Biostrings]{matchPWM}}
#'
#' @export
#' @docType methods
#' @rdname findMotifHits-methods
#'
#' @examples
#' \dontrun{
#' findMotifHits(pwm, "sequences.fa")
#' }
setGeneric(name = "findMotifHits",
           def = function(query, subject, ...) {
               standardGeneric("findMotifHits")
           },
           valueClass = "GRanges")

#' @rdname findMotifHits-methods
#' @aliases findMotifHits,PWMatrix,ANY-method
setMethod("findMotifHits",
          c("PWMatrix","ANY"),
          function(query, subject, ...) {
              findMotifHits(TFBSTools::PWMatrixList(query), subject, ...)
          })

#' @rdname findMotifHits-methods
#' @aliases findMotifHits,PWMatrixList,ANY-method
setMethod("findMotifHits",
          c("PWMatrixList","ANY"),
          function(query, subject, ...) {
              # write motifs to file for homer2
              tmpf <- tempfile(fileext = ".motif")
              .dumpPWMsToHomer2File(pwmL = query, fname = tmpf)
              res <- findMotifHits(tmpf, subject, ...)
              unlink(tmpf)
              return(res)
          })

#' @rdname findMotifHits-methods
#' @aliases findMotifHits,ANY,DNAString-method
setMethod("findMotifHits",
          c("ANY","DNAString"),
          function(query, subject, ...) {
              findMotifHits(query, Biostrings::DNAStringSet(subject), ...)
          })

#' @rdname findMotifHits-methods
#' @aliases findMotifHits,ANY,DNAStringSet-method
setMethod("findMotifHits",
          c("ANY","DNAStringSet"),
          function(query, subject, ...) {
              # write sequences to file for homer2
              tmpf <- tempfile(fileext = ".fa")
              Biostrings::writeXStringSet(x = subject, filepath = tmpf,
                                          append = FALSE, compress = FALSE, format = "fastq")
              res <- findMotifHits(query, tmpf, ...)
              unlink(tmpf)
              return(res)
          })

#' @rdname findMotifHits-methods
#' @aliases findMotifHits,character,character-method
setMethod("findMotifHits",
          c("character","character"),
          function(query, subject, min.score, method = c("homer2", "matchPWM"),
                   homerfile = findHomer("homer2"), Ncpu = 1L, ...) {
    stopifnot(is.character(homerfile) && length(homerfile) == 1L && file.exists(homerfile))
    method <- match.arg(method)
    res <- NULL
    if (method == "matchPWM") {
        stop("method='matchPWM' is not implemented yet")

    } else if (method == "homer2") {

        # make sure motifs are in file
        if (!length(query) == 1L || !file.exists(query)) { # motifs are in file
            stop("'query' must be either a character(1), a PWMatrix or a PWMatrixList object")
        }

        # make sure sequences are in file
        if (length(subject) != 1L || !file.exists(subject)) {
            stop("'subject' must be a character(1) pointing to a FASTA sequence file, ",
                 "or a DNAString, DNAStringSet or GRanges object")
        }

        # run homer2
        cmdargs <- sprintf("find -i %s -m %s -offset 1 -strand both -p %d", subject, query, Ncpu)
        res <- system2(command = homerfile, args = cmdargs, stdout = TRUE, stderr = "", wait = TRUE)
        con <- textConnection(res)
        resparsed <- scan(file = con, what = list(seqnames="", start=1L, matchedSeq="", pwmname="", strand="", score=1.0),
                          sep = "\t", quiet = TRUE)
        close(con)
        res <- GenomicRanges::GRanges(seqnames = resparsed$seqnames,
                                      ranges = IRanges(start=resparsed$start, width=nchar(resparsed$matchedSeq)),
                                      strand = resparsed$strand, matchedSeq = resparsed$matchedSeq,
                                      pwmname = resparsed$pwmname, score = resparsed$score)
    }
    return(res)
})
