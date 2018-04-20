#' @importFrom Biostrings matchPWM DNA_BASES DNAStringSet
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom TFBSTools Matrix PWMatrix PWMatrixList name
NULL

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
#'               file with DNA sequences in FASTA format, or with the sequence(s)})}
#'         \item{\code{DNAString}}{ with a single sequence}
#'         \item{\code{DNAStringSet}}{ with several sequences}
#'         \item{\code{GRanges}}{ object with the genomic coordinates
#'             of the sequences to be searched (.}
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
#' pwm <- TFBSTools::PWMatrix(ID="mypwm", profileMatrix=matrix(rep(c(1,-1,-1,-1),3), nrow=4, dimnames=list(c("A","C","G","T"))))
#' findMotifHits(pwm, "CCCCCAAACCCCC")
setGeneric(name = "findMotifHits",
           def = function(query, subject, min.score, method = c("homer2", "matchPWM"),
                          homerfile = findHomer("homer2"), Ncpu = 1L, ...) {
               method <- match.arg(method)
               if (method == "matchPWM")
                   stop("method='matchPWM' is not implemented yet")
               standardGeneric("findMotifHits")
           },
           valueClass = "GRanges")


#' @rdname findMotifHits-methods
#' @aliases findMotifHits,PWMatrix,character-method
setMethod("findMotifHits",
          c("PWMatrix","character"),
          function(query, subject, ...) {
              findMotifHits(TFBSTools::PWMatrixList(query), subject, ...)
          })

#' @rdname findMotifHits-methods
#' @aliases findMotifHits,PWMatrixList,character-method
setMethod("findMotifHits",
          c("PWMatrix","character"),
          function(query, subject, ...) {
              tmpf <- tempfile(fileext = ".motif")
              fh <- file(tmpf, "wb")
              for (i in seq_along(query)) {
                  pwm <- TFBSTools::Matrix(query[[i]])
                  cat(sprintf(">%s\t%s\t%.2f\n",
                              paste(apply(pwm, 2, function(x) { rownames(pwm)[which.max(x)] }), collapse = ""),
                              TFBSTools::name(query[[i]]), log(2**10)),  file = fh, append = TRUE)
                  write.table(file = fh, t(pwm), row.names = FALSE, col.names = FALSE,
                              sep = "\t", quote = FALSE, append = TRUE)
              }
              close(fh)
              res <- findMotifHits(tmpf, subject, ...)
              unlink(tmpf)
              return(res)
          })

#' @rdname findMotifHits-methods
#' @aliases findMotifHits,character,character-method
setMethod("findMotifHits",
          c("character","character"),
          function(query, subject, min.score, method, homerfile=findHomer("homer2"), Ncpu=1L, ...) {
    stopifnot(is.character(homerfile) && length(homerfile) == 1L && file.exists(homerfile))

    # make sure motifs are in file
    if (length(query) == 1L) {
        if (file.exists(query)) { # motifs are in file
            if (method != "homer2") {
                warning("changed 'method' to 'homer2' for a motif-file query")
                method <- "homer2"
            }
        } else {
            stop("'query' motif file does not exist")
        }
    } else {
        stop("'query' must be either a character(1) or a PWMatrix(List) object")
    }

    # make sure sequences are in file
    if (length(subject) == 1L) { # sequences are in a FASTA file
        if (!file.exists(subject) && sum(strsplit(subject, split = "")[[1]] %in% Biostrings::DNA_BASES) / nchar(subject) >= 0.8) {
            tf <- tempfile(fileext = ".fa")
            on.exit(unlink(subject))
            Biostrings::writeXStringSet(x = DNAStringSet(subject), filepath = tf)
            subject <- tf
        } else {
            stop("'subject' is neither an existing sequence file nor at least 80% DNA_BASES")
        }
    } else {
        tf <- tempfile(fileext = ".fa")
        on.exit(unlink(subject))
        Biostrings::writeXStringSet(x = DNAStringSet(subject), filepath = tf)
        subject <- tf
    }

    # run homer2
    cmdargs <- sprintf("find -i %s -m %s -offset 1 -strand both -p %d", subject, query, Ncpu)
    res <- system2(command = homerfile, args = cmdargs, stdout = TRUE, stderr = "", wait = TRUE)
    con <- textConnection(res)
    resparsed <- scan(file = con, what = list(seqnames="", start=1L, matchedSeq="", pwmname="", strand="", score=1.0),
                      sep = "\t", quiet = TRUE)
    close(con)
    GenomicRanges::GRanges(resparsed$seqnames, IRanges(start=resparsed$start, width=nchar(resparsed$matchedSeq)),
                           strand = resparsed$strand, matchedSeq = resparsed$matchedSeq, pwmname = resparsed$pwmname,
                           score = resparsed$score)
})
