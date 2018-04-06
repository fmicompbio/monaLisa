#' @import GenomicRanges
#' @importFrom tools file_path_as_absolute
#' @importFrom utils write.table
NULL


#' @title Find HOMER script file.
#'
#' @description Find absolute path to HOMER script file.
#'
#' @param scriptfile Name of the script file to search.
#' @param dirs Directory names to look for \code{scriptfile}. If \code{dirs=NULL},
#'     all directories listed in the \code{PATH} environment variable will be
#'     searched.
#'
#' @return Absolute path to \code{scriptfile}, or \code{NA} if none or several were found.
#'
#' @export
findHomerPath <- function(scriptfile = "findMotifsGenome.pl", dirs = NULL) {
    if (is.null(dirs))
        dirs <- strsplit(x = Sys.getenv("PATH"), split = ":")[[1]]
    res <- list.files(path = dirs, pattern = paste0("^",scriptfile,"$"), full.names = TRUE)
    if (length(res) == 1) {
        return(tools::file_path_as_absolute(res))
    } else {
        return(NA)
    }
}


#' @title Dump Jaspar motifs into a HOMER motif file.
#'
#' @description Get PWMs from the \code{JASPAR2016} or \code{JASPAR2018} package
#'     and write them into a HOMER-compatible motif file.
#'
#' @param filename Name of the output file to be created.
#' @param pkg Name of the Jaspar package to use (default: \code{JASPAR2018}).
#' @param opts a search options list used in \code{getMatrixSet}.
#'
#' @return \code{TRUE} if successful.
#'
#' @seealso \code{\link[TFBSTools]{getMatrixSet}} for details on the argument \code{opts}.
#'
#' @export
dumpJaspar <- function(filename, pkg = "JASPAR2018", opts = list(tax_group = "vertebrates")) {
    requireNamespace(pkg)
    requireNamespace("TFBSTools")

    # load PFMs and convert to PWMs
    mdb <- get(pkg)
    siteList <- TFBSTools::getMatrixSet(mdb, opts)
    message("extracted ",length(siteList)," motifs from ",pkg)

    # remark: pseudocount of 4 corresponds to adding 1 to each base count
    #         the following are identical:
    #wm1 <- Matrix(toPWM(siteList[[1]], type="prob", pseudocounts=4.0))
    #wm2 <- Matrix(siteList[[1]])+1 ; wm2 <- t(t(wm2) /colSums(wm2))
    #identical(wm1, wm2)

    message("converting to HOMER format...", appendLF = FALSE)
    fh <- file(filename, "wb")
    for (i in 1:length(siteList)) {
        pwm <- TFBSTools::Matrix(siteList[[i]]) + 1
        pwm <- t(t(pwm) / colSums(pwm))
        tmp.rn <- rownames(pwm)
        pwm <- apply(pwm, 2, function(x){sprintf("%.3f", x)})
        rownames(pwm) <- tmp.rn
        wm.name <- paste(c(TFBSTools::name(siteList[[i]]),
                           paste(TFBSTools::tags(siteList[[i]])$acc, collapse = "::"),
                           TFBSTools::tags(siteList[[i]])$type), collapse = "|")
        cat(sprintf(">%s\t%s\t%.2f\n",
                    paste(apply(pwm, 2, function(x) { rownames(pwm)[which.max(x)] }), collapse = ""),
                    wm.name, log(2**10)),  file = fh, append = TRUE)
        write.table(file = fh, t(pwm), row.names = FALSE, col.names = FALSE,
                    sep = "\t", quote = FALSE, append = TRUE)
        flush(fh)
    }
    close(fh)
    message("done")

    return(TRUE)
}


#' @title Prepare input files for HOMER motif enrichment analysis.
#'
#' @description For each bin, write genomic coordinates for foreground and
#'     background regions into files for HOMER motif enrichment analysis.
#'
#' @param gr A \code{GRanges} object (or an object that can be coerced to one)
#'     with the genomic regions to analyze.
#' @param b A vector of the same length as \code{gr} that groups its elements
#'     into bins (typically a factor).
#' @param outdir A path specifying the folder into which the output files (two
#'     files per unique value of \code{b}) will be written.
#' @param motifFile A file with HOMER formatted PWMs to be used in the enrichment analysis.
#' @param HOMER.path Path and file name of the \code{findMotifsGenome.pl} HOMER script.
#'
#' @details For each bin (unique value of \code{b}) this functions creates two files
#'     in \code{outdir} (\code{outdir/bin_N_foreground.tab} and \code{outdir/bin_N_background.tab},
#'     where \code{N} is the number of the bin and foreground/background correspond to
#'     the ranges that are/are not within the current bin). The files are in the
#'     HOMER peak file format (see http://homer.ucsd.edu/homer/ngs/peakMotifs.html for details).
#'
#'     In addition, a shell script file is created containing the shell commands
#'     to run the HOMER motif enrichment analysis.
#'
#' @return The path and name of the script file to run the HOMER motif enrichment
#'     analysis.
#'
#' @export
prepareHOMER <- function(gr, b, outdir, motifFile, HOMER.path = findHomerPath()) {
    stopifnot(inherits(gr, "GRanges"))
    stopifnot(is.vector(b) && length(b) == length(gr))
    stopifnot(is.character(outdir))
    stopifnot(file.exists(motifFile))
    stopifnot(file.exists(HOMER.path))

    if (file.exists(outdir))
        stop(outdir," already exists - will not overwrite existing folder")

    dir.create(outdir)
}
