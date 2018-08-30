#' @importFrom Biostrings matchPWM DNA_BASES DNAStringSet
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom TFBSTools Matrix PWMatrix PWMatrixList name
#' @importFrom methods .valueClassTest
#' @importFrom utils strcapture
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

# internal function:  read motifs from file in homer2 format into PWMatrixList
.readPWMsFromHomer2File <- function(fname) {
    stopifnot(is.character(fname) && length(fname) == 1L && file.exists(fname))
    lns <- readLines(fname)
    i <- grep("^>", lns)
    df <- strcapture(pattern = "^>([^\t]+)\t([^\t]+)\t([-0-9.]+).*$", x = lns[i],
                     proto = data.frame(id = character(), name = character(), cutoffScore = numeric()))
    df$id <- as.character(df$id)
    df$name <- as.character(df$name)
    w <- diff(c(i, length(lns) + 1L)) - 1L
    mL <- split(lns[-i], rep(seq_along(i), w))
    names(mL) <- df$name
    mL <- lapply(mL, function(x) {
        m <- do.call(cbind, lapply(x, function(xx) as.numeric(unlist(strsplit(xx, "\t", fixed = TRUE)))))
        rownames(m) <- c("A", "C", "G", "T")
        m
    })
    pwmL <- lapply(seq_along(mL), function(j) PWMatrix(ID = df[j,"id"], name = df[j,"name"], profileMatrix = mL[[j]]))
    pwmL <- do.call(TFBSTools::PWMatrixList, pwmL)
    names(pwmL) <- df$name
    return(pwmL)
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

# method hierarchy:
# 1. findMotifHits,character,character        --> if method=="homer2", use homer2; else --> 9.
# 2. findMotifHits,character,DNAString        --> 3.
# 3. findMotifHits,character,DNAStringSet     --> if method=="homer2" --> 1.     ; else --> 9.
# 4. findMotifHits,PWMatrix,character         --> 7.
# 5. findMotifHits,PWMatrix,DNAString         --> 9.
# 6. findMotifHits,PWMatrix,DNAStringSet      --> 9.
# 7. findMotifHits,PWMatrixList,character     --> if method=="homer2" --> 1.     ; else --> 9.
# 8. findMotifHits,PWMatrixList,DNAString     --> 9.
# 9. findMotifHits,PWMatrixList,DNAStringSet


# 1.
#' @rdname findMotifHits-methods
#' @aliases findMotifHits,character,character-method
setMethod("findMotifHits",
          c("character","character"),
          function(query, subject, min.score, method = c("homer2", "matchPWM"),
                   homerfile = findHomer("homer2"), Ncpu = 1L, ...) {
              method <- match.arg(method)
              stopifnot(method == "matchPWM" || (method == "homer2" && is.character(homerfile) && length(homerfile) == 1L && file.exists(homerfile)))
              if (method == "matchPWM") {
                  # read in motifs from query into PMatrixList
                  pwmL <- .readPWMsFromHomer2File(query)
                  # read in sequences from subject into DNAStringSet
                  seqs <- Biostrings::readDNAStringSet(filepath = subject)
                  # call next method
                  findMotifHits(query = pwmL, subject = seqs, min.score = min.score, method = method, Ncpu = Ncpu, ...)

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
                  resparsed <- scan(file = con, what = list(seqnames = "", start = 1L, matchedSeq = "",
                                                            pwmname = "", strand = "", score = 1.0),
                                    sep = "\t", quiet = TRUE)
                  close(con)
                  sl <- Biostrings::fasta.seqlengths(subject)
                  gr <- GenomicRanges::GRanges(seqnames = resparsed$seqnames,
                                               ranges = IRanges(start = resparsed$start, width = nchar(resparsed$matchedSeq)),
                                               strand = resparsed$strand, matchedSeq = resparsed$matchedSeq,
                                               pwmname = resparsed$pwmname, score = resparsed$score,
                                               seqlengths = sl)
                  sort(gr)
              } else {
                  stop("unknown 'method': ", method)
              }
          })

# 2.
#' @rdname findMotifHits-methods
#' @aliases findMotifHits,character,DNAString-method
setMethod("findMotifHits",
          c("character","DNAString"),
          function(query, subject, ...) {
              findMotifHits(query, Biostrings::DNAStringSet(subject), ...)
          })

# 3.
#' @rdname findMotifHits-methods
#' @aliases findMotifHits,character,DNAStringSet-method
setMethod("findMotifHits",
          c("character","DNAStringSet"),
          function(query, subject, min.score, method = c("homer2", "matchPWM"),
                   homerfile = findHomer("homer2"), Ncpu = 1L, ...) {
              method <- match.arg(method)
              if (method == "matchPWM") {
                  # read in motifs from query into PMatrixList
                  pwmL <- .readPWMsFromHomer2File(query)
                  # call next method
                  findMotifHits(query = pwmL, subject = subject, min.score = min.score,
                                method = method, homerfile = homerfile, Ncpu = Ncpu, ...)

              } else if (method == "homer2") {
                  # write sequences to file
                  tmpf <- tempfile(fileext = ".fa")
                  Biostrings::writeXStringSet(x = subject, filepath = tmpf,
                                              append = FALSE, compress = FALSE, format = "fasta")
                  res <- findMotifHits(query = query, subject = tmpf, min.score = min.score,
                                       method = method, homerfile = homerfile, Ncpu = Ncpu, ...)
                  unlink(tmpf)
                  return(res)
              } else {
                  stop("unknown 'method': ", method)
              }
          })

# 4.
#' @rdname findMotifHits-methods
#' @aliases findMotifHits,PWMatrix,character-method
setMethod("findMotifHits",
          c("PWMatrix","character"),
          function(query, subject, ...) {
              findMotifHits(TFBSTools::PWMatrixList(query), subject, ...)
          })

# 5.
#' @rdname findMotifHits-methods
#' @aliases findMotifHits,PWMatrix,DNAString-method
setMethod("findMotifHits",
          c("PWMatrix","DNAString"),
          function(query, subject, ...) {
              findMotifHits(TFBSTools::PWMatrixList(query), Biostrings::DNAStringSet(subject), ...)
          })

# 6.
#' @rdname findMotifHits-methods
#' @aliases findMotifHits,PWMatrix,DNAStringSet-method
setMethod("findMotifHits",
          c("PWMatrix","DNAStringSet"),
          function(query, subject, ...) {
              findMotifHits(TFBSTools::PWMatrixList(query), subject, ...)
          })

# 7.
#' @rdname findMotifHits-methods
#' @aliases findMotifHits,PWMatrixList,character-method
setMethod("findMotifHits",
          c("PWMatrixList","character"),
          function(query, subject, min.score, method = c("homer2", "matchPWM"),
                   homerfile = findHomer("homer2"), Ncpu = 1L, ...) {
              method <- match.arg(method)
              if (method == "matchPWM") {
                  # read in sequences from subject into DNAStringSet
                  seqs <- Biostrings::readDNAStringSet(filepath = subject)
                  # call next method
                  findMotifHits(query = query, subject = seqs, min.score = min.score,
                                method = method, homerfile = homerfile, Ncpu = Ncpu, ...)

              } else if (method == "homer2") {
                  # write motifs to file
                  tmpf <- tempfile(fileext = ".motif")
                  if (is.character(min.score) && grepl(pattern = "^[0-9.]+[%]$", x = min.score))
                      .dumpPWMsToHomer2File(pwmL = query, fname = tmpf, relscore = as.numeric(sub("[%]$","",min.score))/100)
                  else if (is.numeric(min.score))
                      .dumpPWMsToHomer2File(pwmL = query, fname = tmpf, absscore = min.score)
                  else
                      stop("wrong type of 'min.score': ",min.score)
                  res <- findMotifHits(query = tmpf, subject = subject, min.score = min.score,
                                       method = method, homerfile = homerfile, Ncpu = Ncpu, ...)
                  unlink(tmpf)
                  return(res)
              } else {
                  stop("unknown 'method': ", method)
              }
          })

# 8.
#' @rdname findMotifHits-methods
#' @aliases findMotifHits,PWMatrixList,DNAString-method
setMethod("findMotifHits",
          c("PWMatrixList","DNAString"),
          function(query, subject, ...) {
              findMotifHits(query, Biostrings::DNAStringSet(subject), ...)
          })

# 9.
#' @rdname findMotifHits-methods
#' @aliases findMotifHits,PWMatrixList,DNAStringSet-method
setMethod("findMotifHits",
          c("PWMatrixList","DNAStringSet"),
          function(query, subject, min.score, method = c("homer2", "matchPWM"),
                   homerfile = findHomer("homer2"), Ncpu = 1L, ...) {
              method <- match.arg(method)
              if (method == "matchPWM") {
                  # search motifs using matchPWM (TFBSTools::searchSeq uses matchPWM internally, but has a nicer interface)
                  tmp <- TFBSTools::searchSeq(x = query, subject = subject, min.score = min.score, mc.cores = Ncpu)
                  df <- as(tmp, "DataFrame")
                  sl <- width(subject); names(sl) <- names(subject)
                  gr <- GenomicRanges::GRanges(seqnames = df$seqnames,
                                               ranges = IRanges(start = df$start, end = df$end),
                                               strand = df$strand, matchedSeq = df$siteSeqs,
                                               pwmname = df$TF, score = df$absScore,
                                               seqlengths = sl)
                  sort(gr)

              } else if (method == "homer2") {
                  # write motifs to file
                  tmpf1 <- tempfile(fileext = ".motif")
                  if (is.character(min.score) && grepl(pattern = "^[0-9.]+[%]$", x = min.score))
                      .dumpPWMsToHomer2File(pwmL = query, fname = tmpf1, relscore = as.numeric(sub("[%]$","",min.score))/100)
                  else if (is.numeric(min.score))
                      .dumpPWMsToHomer2File(pwmL = query, fname = tmpf1, absscore = min.score)
                  else
                      stop("wrong type of 'min.score': ",min.score)
                  # write sequences to file
                  tmpf2 <- tempfile(fileext = ".fa")
                  Biostrings::writeXStringSet(x = subject, filepath = tmpf2,
                                              append = FALSE, compress = FALSE, format = "fasta")
                  # call next method
                  res <- findMotifHits(query = tmpf1, subject = tmpf2, min.score = min.score,
                                       method = method, homerfile = homerfile, Ncpu = Ncpu, ...)
                  unlink(c(tmpf1, tmpf2))
                  return(res)
              } else {
                  stop("unknown 'method': ", method)
              }
          })
