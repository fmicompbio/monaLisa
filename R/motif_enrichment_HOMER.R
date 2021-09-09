#' @import GenomicRanges
#' @importFrom tools file_path_as_absolute
#' @importFrom utils read.delim write.table getFromNamespace
#' @importFrom stats p.adjust
#' @importFrom methods as
NULL


#' @title Find HOMER script file.
#'
#' @description Find absolute path to HOMER script file.
#'
#' @param homerfile Name of the script file to search.
#' @param dirs Directory names to look for \code{homerfile}. If \code{dirs=NULL},
#'     all directories listed in the \code{PATH} environment variable will be
#'     searched.
#' 
#' @details In addition to \code{dirs}, \code{findHomer} will also look in the
#'     directory provided in the environment variable \code{MONALISA_HOMER}. 
#'
#' @return Absolute path to \code{homerfile}, or \code{NA} if none or several were found.
#' 
#' @examples 
#' homer_path <- findHomer()
#' 
#' @export
findHomer <- function(homerfile = "findMotifsGenome.pl", dirs = NULL) {
    if (is.null(dirs))
        dirs <- strsplit(x = Sys.getenv("PATH"), split = ":")[[1]]
    if (!is.na(monalisa_homer <- Sys.getenv("MONALISA_HOMER", unset = NA)))
        dirs <- c(strsplit(x = monalisa_homer, split = ":")[[1]], dirs)
    dirs <- unique(dirs)
    res <- list.files(path = dirs, pattern = paste0("^",homerfile,"$"),
                      full.names = TRUE, recursive = FALSE)
    if (length(res) == 1) {
        return(tools::file_path_as_absolute(res))
    } else {
        return(NA)
    }
}


#' @title Dump Jaspar motifs into a HOMER motif file.
#'
#' @description Get motifs from a Jaspar database package (e.g. \code{JASPAR2020})
#'     and write them into a HOMER-compatible motif file as positional probability
#'     matrices.
#'
#' @param filename Name of the output file to be created.
#' @param pkg Name of the Jaspar package to use (default: \code{JASPAR2020}).
#' @param opts A list with search options used in
#'     \code{\link[TFBSTools]{getMatrixSet}}. By default, only vertebrate motifs
#'     are included in the output using
#'     \code{opts = list(tax_group = "vertebrates")}.
#' @param pseudocount A numerical scalar with the pseudocount to be added to
#'     each element of the position frequency matrix extracted from Jaspar,
#'     before its conversion to a position probability matrix (default: 1.0).
#' @param relScoreCutoff Currently ignored. numeric(1) in [0,1] that sets the default motif
#'     log-odds score cutof to relScoreCutoff * maximal score for each PWM
#'     (default: 0.8).
#' @param verbose A logical scalar. If \code{TRUE}, print progress messages.
#'
#' @return \code{TRUE} if successful.
#'
#' @seealso \code{\link[TFBSTools]{getMatrixSet}} for details on the argument \code{opts}.
#'     \code{\link{homerToPFMatrixList}} to read a file with HOMER-formatted
#'     motifs into a \code{\link[TFBSTools]{PFMatrixList}}.
#'
#' @examples 
#' dumpJaspar(filename = tempfile(), pkg = "JASPAR2020", 
#'            opts = list(ID = c("MA0006.1")))
#' 
#' @importFrom utils getFromNamespace
#' @importFrom TFBSTools getMatrixSet Matrix ID name
#'
#' @export
dumpJaspar <- function(filename,
                       pkg = "JASPAR2020",
                       opts = list(tax_group = "vertebrates"),
                       pseudocount = 1,
                       relScoreCutoff = 0.8,
                       verbose = FALSE) {
    .assertScalar(x = filename, type = "character")
    stopifnot(!file.exists(filename))
    .assertScalar(x = pseudocount, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = relScoreCutoff, type = "numeric", rngIncl = c(0, 1))
    if ("matrixtype" %in% names(opts) && opts[["matrixtype"]] != "PFM") {
      stop("opts[['matrixtype']] must be set to 'PFM'")
    } else {
      opts[["matrixtype"]] <- "PFM"
    }
    .assertScalar(x = verbose, type = "logical")

    # load PFMs and convert to HOMER format (base probabilties)
    requireNamespace(pkg)
    mdb <- utils::getFromNamespace(pkg, ns = pkg)
    siteList <- TFBSTools::getMatrixSet(mdb, opts)
    if (verbose)
        message("extracted ",length(siteList)," motifs from ",pkg)

    # remark: pseudocount of 4 corresponds to adding 1 to each base count
    #         the following are identical:
    #wm1 <- Matrix(toPWM(siteList[[1]], type="prob", pseudocounts=4.0))
    #wm2 <- Matrix(siteList[[1]])+1 ; wm2 <- t(t(wm2) /colSums(wm2))
    #identical(wm1, wm2)

    if (verbose)
        message("converting to HOMER format...", appendLF = FALSE)
    fh <- file(filename, "wb")
    for (i in seq_len(length(siteList))) {
        ppm <- TFBSTools::Matrix(siteList[[i]]) + pseudocount
        ppm <- t(t(ppm) / colSums(ppm))
        tmp.rn <- rownames(ppm)
        #scorecut <- relScoreCutoff * sum(log(apply(ppm, 2, max) / 0.25))
        scorecut <- log(2**10) # use constant cutoff for all motifs (ignore relScoreCutoff)
        ppm <- apply(ppm, 2, function(x){sprintf("%.3f", x)})
        rownames(ppm) <- tmp.rn
        wm.name <- paste0(TFBSTools::ID(siteList[[i]]), ":::", TFBSTools::name(siteList[[i]]))

        # the -10 is added so that the motif file has 4 columns,
        # which is needed to run compareMotifs.pl for weight matrix clustering
        # (4th column not used, bug in compareMotifs.pl, I think)
        cat(sprintf(">%s\t%s\t%.2f\t-10\n",
                    paste(apply(ppm, 2, function(x) {
                        rownames(ppm)[which.max(x)]
                    }), collapse = ""),
                    wm.name, scorecut),  file = fh, append = TRUE)
        write.table(file = fh, t(ppm), row.names = FALSE, col.names = FALSE,
                    sep = "\t", quote = FALSE, append = TRUE)
        flush(fh)
    }
    close(fh)
    if (verbose)
        message("done")

    return(TRUE)
}


#' @title Read a HOMER motif file and create a PFMatrixList
#'
#' @description Read motifs from a file in HOMER format and create
#'     a PFMatrixList from them.
#'
#' @param filename Name of the input file with HOMER-formatted motifs.
#' @param n The number of observations (multiplied with base frequencies to
#'     create the number of observed bases at each position).
#'
#' @return A \code{\link[TFBSTools]{PFMatrixList}} with motifs from the file.
#'
#' @examples 
#' library(JASPAR2020)
#' optsL <- list(ID = c("MA0006.1"))
#' pfm1 <- TFBSTools::getMatrixSet(JASPAR2020, opts = optsL)
#' TFBSTools::Matrix(pfm1)
#' 
#' tmpfn <- tempfile()
#' dumpJaspar(filename = tmpfn, pkg = "JASPAR2020", opts = optsL)
#' pfm2 <- homerToPFMatrixList(tmpfn)
#' TFBSTools::Matrix(pfm2)
#' 
#' unlink(tmpfn)
#' 
#' @seealso \code{\link{dumpJaspar}} for writing motifs from a Jaspar database
#'     package into a file in HOMER format.
#'
#' @importFrom utils getFromNamespace
#' @importFrom TFBSTools PFMatrix PFMatrixList
#'
#' @export
homerToPFMatrixList <- function(filename, n = 100L) {
    .assertScalar(x = filename, type = "character")
    stopifnot(file.exists(filename))
    .assertScalar(x = n, type = "numeric", rngExcl = c(0, Inf))

    # parse HOMER motif file
    tmp <- readLines(filename)
    g <- grep(">", tmp)
    tmp[g] <- sub("^>","",tmp[g])
    fields <- strsplit(tmp[g], "\t", fixed = TRUE)
    L <- lapply(seq_along(g), function(i) {
        cons <- fields[[i]][1]
        nm <- fields[[i]][2]
        log2cut <- round(as.numeric(fields[[i]][3]) / log(2), 2)
        s <- g[i] + 1L
        e <- if (i == length(g)) length(tmp) else g[i + 1L] - 1L
        m <- round(n * do.call(cbind, lapply(strsplit(tmp[s:e], "\t", fixed = TRUE), as.numeric)), 0)
        rownames(m) <- c("A", "C", "G", "T")
        TFBSTools::PFMatrix(ID = nm, name = nm, profileMatrix = m,
                            tags = list(log2cut = log2cut, comment = "imported from HOMER motif file"))
    })

    do.call(PFMatrixList, c(L, list(use.names = TRUE)))
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
#' @param genomedir Directory containing sequence files in Fasta format (one per chromosome).
#' @param outdir A path specifying the folder into which the output files (two
#'     files per unique value of \code{b}) will be written.
#' @param motifFile A file with HOMER formatted PWMs to be used in the enrichment analysis.
#' @param homerfile Path and file name of the \code{findMotifsGenome.pl} HOMER script.
#' @param regionsize The peak size to use in HOMER (\code{"given"} keeps the coordinate
#'     region, an integer value will keep only that many bases in the region center).
#' @param Ncpu Number of parallel threads that HOMER can use.
#' @param verbose A logical scalar. If \code{TRUE}, print progress messages.
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
#' @examples
#' # prepare genome directory (here: one dummy chromosome)
#' genomedir <- tempfile()
#' dir.create(genomedir)
#' writeLines(c(">chr1", "ATGCATGCATCGATCGATCGATCGTACGTA"),
#'            file.path(genomedir, "chr1.fa"))
#' 
#' # prepare motif file, regions and bins
#' motiffile <- tempfile()
#' dumpJaspar(filename = motiffile, pkg = "JASPAR2020", 
#'            opts = list(ID = c("MA0006.1")))
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1:4, width = 4))
#' b <- bin(1:4, nElements = 2)
#'
#' # create dummy file (should point to local Homer installation)
#' homerfile <- file.path(tempdir(), "findMotifsGenome.pl")
#' writeLines("dummy", homerfile)
#' 
#' # run prepareHomer
#' outdir <- tempfile()
#' prepareHomer(gr = gr, b = b, genomedir = genomedir,
#'              outdir = outdir, motifFile = motiffile,
#'              homerfile = homerfile, verbose = TRUE) 
#' list.files(outdir)
#' 
#' # clean up example
#' unlink(c(genomedir, motiffile, homerfile, outdir))
#'
#' @importFrom GenomicRanges start end strand
#' @export
prepareHomer <- function(gr, b, genomedir, outdir, motifFile,
                         homerfile = findHomer(), regionsize = "given",
                         Ncpu = 2L, verbose = FALSE) {
    if (!inherits(gr, "GRanges"))
        as(gr, "GRanges")
    if (!is.factor(b))
        b <- factor(b, levels = unique(b))
    .assertVector(x = b, type = "factor", len = length(gr))
    .assertScalar(x = genomedir, type = "character")
    .assertScalar(x = outdir, type = "character")
    .assertScalar(x = motifFile, type = "character")
    .assertScalar(x = homerfile, type = "character")
    stopifnot(exprs = {
        file.exists(genomedir)
        file.exists(motifFile)
        file.exists(homerfile)
        regionsize == "given" || (is.numeric(regionsize) && length(regionsize) == 1L && regionsize > 0)
    })
    .assertScalar(x = verbose, type = "logical")

    if (file.exists(outdir))
        stop(outdir," already exists - will not overwrite existing folder")
    dir.create(outdir)

    homerFile <- file.path(outdir, "run.sh")
    fh <- file(homerFile, "w")

    if (verbose)
        message("creating foreground/background region files for HOMER")
    for (i in seq_len(nlevels(b))) {
        bn <- levels(b)[i]
        if (verbose)
            message("  bin ",bn)

        fgfile  <- sprintf("%s/bin_%03d_foreground.tab", outdir, i)
        bgfile  <- sprintf("%s/bin_%03d_background.tab", outdir, i)
        outputf <- sprintf("%s/bin_%03d_output", outdir, i)

        tmp.gr <- gr[b == bn]
        write.table(file = fgfile,
                    data.frame(seq_along(tmp.gr), as.character(seqnames(tmp.gr)),
                               start(tmp.gr) - 1, end(tmp.gr),
                               ifelse(as.character(strand(tmp.gr)) == "-", "-", "+")),
                    sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

        tmp.gr <- gr[b != bn]
        write.table(file = bgfile,
                    data.frame(seq_along(tmp.gr), as.character(seqnames(tmp.gr)),
                               start(tmp.gr) - 1, end(tmp.gr), 
                               ifelse(as.character(strand(tmp.gr)) == "-", "-", "+")),
                    sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

        cat(sprintf("%s %s %s %s -bg %s -nomotif -p %d -size %s -mknown %s\n",
                    homerfile, fgfile, genomedir, outputf, bgfile, Ncpu,
                    as.character(regionsize), motifFile), file = fh, append = TRUE)
    }
    close(fh)

    return(homerFile)
}

#' @title load output from HOMER findMotifsGenome.pl into R
#'
#' @description Parse HOMER output files into R data structures.
#'
#' @param infiles HOMER output files to be parsed.
#' @param pseudocount.log2enr A numerical scalar with the pseudocount to add to
#'   foreground and background counts when calculating log2 motif enrichments
#' @param pseudofreq.pearsonResid A numerical scalar with the pseudo-frequency
#'   to add to background frequencies when calculating Pearson residuals.
#'   The value needs to be in [0,1] and corresponds to the minimal expected
#'   fraction of background sequences that contain at least one motif hit.
#' @param p.adjust.method A character scalar selecting the p value adjustment
#'   method (used in \code{\link[stats]{p.adjust}}).
#'
#' @return A list of four components (\code{negLog10P}, \code{negLog10Padj},
#'     \code{pearsonResid} and \code{log2enr}), containing each a motif (rows)
#'     by bin (columns) matrix with raw -log10 P values, -log10 adjusted P values,
#'     and motif enrichments as Pearson residuals (\code{pearsonResid}) and as
#'     log2 ratios (\code{log2enr}).
#'
#' @examples 
#' outfile <- system.file("extdata", "homer_output.txt.gz", package = "monaLisa")
#' res <- parseHomerOutput(infiles = c(bin1 = outfile))
#' head(res$negLog10P)
#' 
#' @importFrom stats p.adjust
#'
#' @export
parseHomerOutput <- function(infiles,
                             pseudocount.log2enr = 8,
                             pseudofreq.pearsonResid = 0.001,
                             p.adjust.method = "BH") {
    stopifnot(all(file.exists(infiles)))
    .assertScalar(x = pseudocount.log2enr, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = pseudofreq.pearsonResid, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = p.adjust.method, type = "character", validValues = stats::p.adjust.methods)

    tabL <- lapply(infiles, read.delim)
    mnms <- sort(tabL[[1]][, 1])
    names(tabL) <- if (is.null(names(infiles))) infiles else names(infiles)
    P <- do.call(cbind, lapply(tabL, function(x) -x[match(mnms, x[, 1]), 4] / log(10)))
    # ... Pearson residuals
    #     assuming expTF to be a Binomial random variable, with
    #       mean     = N_fg * p_bg
    #       variance = N_fg * p_bg * (1 - p_bg)
    #     idea: each sequence is one trial with outcomes "hit" or "no hit"
    #     (zoops), with a constant "hit" rate within a sequence set
    #     (p_bg in the background); expTF thus follows a binomial distribution
    #     with the number of trials corresponding to the (weighted) number of
    #     sequences (e.g. N_fg)
    presid <- do.call(cbind, lapply(tabL, function(tab) {
        nTotFg <- as.numeric(gsub("\\S+\\.(\\d+)\\.", "\\1", colnames(tab)[6])) #total number of target (foreground) sequences
        fracBgWithMotif <-
          pmin(1, as.numeric(sub("%$","", tab[, 9])) / 100 +
                 pseudofreq.pearsonResid)
        obsTF <- tab[, 6]
        expTF <- nTotFg * fracBgWithMotif
        enr <- (obsTF - expTF) / sqrt(expTF * (1 - fracBgWithMotif))
        enr[ is.na(enr) ] <- 0
        enr[match(mnms, tab[, 1])]
    }))
    log2enr <- do.call(cbind, lapply(tabL, function(tab){
        numFgBgWithHits <- tab[, c(6, 8)] # number of target seqs and bg seqs with motif
        nTot <- as.numeric(gsub("\\S+\\.(\\d+)\\.", "\\1", colnames(numFgBgWithHits))) #total number of target and background sequences
        numFgBgWithHitsNorm <- t(min(nTot) * t(numFgBgWithHits) / nTot) # scale to smaller number (usually number of target sequences)
        numFgBgWithHitsNormLog <- log2(numFgBgWithHitsNorm + pseudocount.log2enr)
        lenr <- numFgBgWithHitsNormLog[, 1] - numFgBgWithHitsNormLog[, 2]
        lenr[match(mnms, tab[, 1])]
    }))

    sumFgWgt <- do.call(cbind, lapply(tabL, function(tab) tab[match(mnms, tab[, 1]), 6]))
    sumBgWgt <- do.call(cbind, lapply(tabL, function(tab) tab[match(mnms, tab[, 1]), 8]))

    totWgt <- do.call(rbind, lapply(tabL, function(tab) {
        numFgBgWithHits <- tab[, c(6, 8)] # number of target seqs and bg seqs with motif
        nTot <- as.numeric(gsub("\\S+\\.(\\d+)\\.", "\\1", colnames(numFgBgWithHits))) #total number of target and background sequences
        nTot
    }))

    padj <- matrix(-log10(p.adjust(as.vector(10^(-P)), method = p.adjust.method)),
                  nrow = nrow(P), dimnames = dimnames(P))
    padj[which(padj == Inf, arr.ind = TRUE)] <- max(padj[is.finite(padj)])
    
    rownames(P) <- rownames(padj) <- rownames(presid) <- rownames(log2enr) <-
        rownames(sumFgWgt) <- rownames(sumBgWgt) <- mnms
    
    return(list(negLog10P = P, negLog10Padj = padj,
                pearsonResid = presid, log2enr = log2enr,
                sumForegroundWgtWithHits = sumFgWgt, 
                sumBackgroundWgtWithHits = sumBgWgt,
                totalWgtForeground = totWgt[, 1],
                totalWgtBackground = totWgt[, 2]))
}

# internal function:  Check if all the HOMER output files already exist and if the run was successful.
# - needs the motifFile and outdir inputs to calcBinnedMotifEnrHomer.
#' @importFrom utils read.table
.checkHomerRun <- function(motifFile, outdir, nbins){

    # Checks
    # --> check resultsfile is a txt file

    # get knownResults.txt files
    out_files <- dir(path = outdir, pattern = "knownResults.txt",
                     full.names = TRUE, recursive = TRUE, ignore.case = FALSE)

    # case 1: out_files is empty, so return FALSE
    if (isEmpty(out_files)) {
      return(FALSE)
    }

    # case 2: at least one file exists, so we check HOMER ran correctly and that the number of output files matches the number of bins

      # get motif names from motifFile
      lns <- readLines(motifFile)
      motifs_motifFile <- unlist(lapply(strsplit(lns[grep("^>", lns)], "\t"), "[", 2))

      # get motif names from resultsfile
      df_list <- lapply(as.list(out_files), function(f){read.table(f, header = FALSE, sep = "\t", skip = 1)}) # skip the column names
      motifs_outdir_list <- lapply(df_list, function(df){ as.character(df[ ,1]) })
      # check the completeness of the run and that the number of output files equals number of bins
      all(vapply(motifs_outdir_list, function(x){all(motifs_motifFile %in% x)}, NA)) & (length(out_files) == nbins)

}

#' @title Prepare and run HOMER motif enrichment analysis.
#'
#' @description Run complete HOMER motif enrichment analysis, consisting of
#'     calls to \code{\link{prepareHomer}}, \code{\link[base]{system2}} and
#'     \code{\link{parseHomerOutput}}.
#'
#' @param gr A \code{GRanges} object (or an object that can be coerced to one)
#'     with the genomic regions to analyze.
#' @param b A vector of the same length as \code{gr} that groups its elements
#'     into bins (typically a factor, such as the one returned by \code{\link{bin}}).
#' @param genomedir Directory containing sequence files in Fasta format (one per chromosome).
#' @param outdir A path specifying the folder into which the output files will be written.
#' @param motifFile A file with HOMER formatted PWMs to be used in the enrichment analysis.
#' @param homerfile Path and file name of the \code{findMotifsGenome.pl} HOMER script.
#' @param regionsize The peak size to use in HOMER (\code{"given"} keeps the coordinate
#'     region, an integer value will keep only that many bases in the region center).
#' @param pseudocount.log2enr A numerical scalar with the pseudocount to add to
#'   foreground and background counts when calculating log2 motif enrichments
#' @param pseudofreq.pearsonResid A numerical scalar with the pseudo-frequency
#'   to add to background frequencies when calculating Pearson residuals.
#'   The value needs to be in [0,1] and corresponds to the minimal expected
#'   fraction of background sequences that contain at least one motif hit.
#' @param p.adjust.method A character scalar selecting the p value adjustment
#'   method (used in \code{\link[stats]{p.adjust}}).
#' @param Ncpu Number of parallel threads that HOMER can use.
#' @param verbose A logical scalar. If \code{TRUE}, print progress messages.
#' @param verbose.Homer A logical scalar. If \code{TRUE}, print the console
#'   output when running Homer.
#'
#' @seealso The functions that are wrapped: \code{\link{prepareHomer}},
#'     \code{\link[base]{system2}} and \code{\link{parseHomerOutput}},
#'     \code{\link{bin}} for binning of regions
#'
#' @return A \code{SummarizedExperiment} object with motifs in rows and bins
#'   in columns, containing six assays: \itemize{
#'   \item{negLog10P}{: -log10 P values}
#'   \item{negLog10Padj}{: -log10 adjusted P values}
#'   \item{pearsonResid}{: motif enrichments as Pearson residuals}
#'   \item{log2enr}{: motif enrichments as log2 ratios}
#'   \item{sumForegroundWgtWithHits}{: Sum of foreground sequence weights
#'     in a bin that have motif hits}
#'   \item{sumBackgroundWgtWithHits}{: Sum of background sequence weights
#'     in a bin that have motif hits}
#' }
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom TFBSTools ID name Matrix
#'
#' @export
calcBinnedMotifEnrHomer <- function(gr, b, genomedir, outdir, motifFile,
                                    homerfile = findHomer(),
                                    regionsize = "given",
                                    pseudocount.log2enr = 8,
                                    pseudofreq.pearsonResid = 0.001,
                                    p.adjust.method = "BH",
                                    Ncpu = 2L,
                                    verbose = FALSE,
                                    verbose.Homer = FALSE) {
    if (!inherits(gr, "GRanges"))
        gr <- as(gr, "GRanges")
    if (!is.factor(b))
        b <- factor(b, levels = unique(b))
    .assertVector(x = b, type = "factor", len = length(gr))
    .assertScalar(x = genomedir, type = "character")
    .assertScalar(x = outdir, type = "character")
    .assertScalar(x = motifFile, type = "character")
    .assertScalar(x = homerfile, type = "character")
    stopifnot(exprs = {
        file.exists(genomedir)
        file.exists(motifFile)
        file.exists(homerfile)
        regionsize == "given" || (is.numeric(regionsize) && length(regionsize) == 1L && regionsize > 0)
    })
    .assertScalar(x = Ncpu, type = "numeric", rngIncl = c(1, Inf))
    .assertScalar(x = verbose, type = "logical")
    .assertScalar(x = verbose.Homer, type = "logical")
    
    ## ... check if the HOMER output is already there for all bins and if it ran completely:
    ## ... ... If yes, go to the 'parse output step', otherwise run homer and check again
    if (.checkHomerRun(motifFile = motifFile, outdir = outdir, nbins = nlevels(b))) {

        if (verbose)
            message("\nHOMER output files already exist, using existing files...")

    } else {

      ## ... case: all/some files exist and/or HOMER didn't run correctly: warn the User to remove all existing files and rerun
      if (any(file.exists(dir(path = outdir, pattern = "knownResults.txt",
                              full.names = TRUE, recursive = TRUE, ignore.case = FALSE)))) {

          stop("\nThere are existing 'knownResults.txt' file(s) in outdir. ",
               "There may be missing 'knownResults.txt' files for some bins ",
               "and/or the existing files are incomplete (cases where the HOMER ",
               "run failed). Please delete these files and rerun 'calcBinnedMotifEnrHomer'.")

      }

      ## ... prepare
      if (verbose)
          message("\npreparing input files...")
      runfile <- prepareHomer(gr = gr, b = b, genomedir = genomedir, outdir = outdir,
                              motifFile = motifFile, homerfile = homerfile,
                              regionsize = regionsize, Ncpu = Ncpu)

      ## ... run
      if (verbose)
          message("\nrunning HOMER...")
      system2(command = "sh", args = runfile,
              stdout = ifelse(verbose.Homer, "", FALSE),
              stderr = ifelse(verbose.Homer, "", FALSE),
              env = paste0("PATH=", dirname(homerfile), ":", Sys.getenv("PATH"), ";"))

      ## ... check HOMER ran correctly
      if (!.checkHomerRun(motifFile = motifFile, outdir = outdir, nbins = length(levels(b)))) {
          stop("HOMER output wasn't complete. Try running again.")
      }

    }

    ## ... parse output
    resfiles <- sprintf("%s/bin_%03d_output/knownResults.txt", outdir, seq_along(levels(b)))
    names(resfiles) <- levels(b)
    resL <- parseHomerOutput(infiles = resfiles,
                             pseudocount.log2enr = pseudocount.log2enr,
                             p.adjust.method = p.adjust.method)

    ## ... create SummarizedExperiment
    ## ... ... reorder parsed output according to motifs in motifFile
    pfms <- homerToPFMatrixList(motifFile)
    morder <- TFBSTools::name(pfms)
    o <- match(morder, rownames(resL[[1]]))
    assayL <- lapply(resL[seq_len(6)], function(x) x[o, ])
    ## ... ... colData
    if (is.null(attr(b, "breaks"))) {
        binL <- binH <- rep(NA, nlevels(b))
    } else {
        binL <- attr(b, "breaks")[-(nlevels(b) + 1)]
        binH <- attr(b, "breaks")[-1]
    }
    cdat <- S4Vectors::DataFrame(bin.names = levels(b),
                                 bin.lower = binL,
                                 bin.upper = binH,
                                 bin.nochange = seq.int(nlevels(b)) %in% getZeroBin(b),
                                 totalWgtForeground = resL$totalWgtForeground,
                                 totalWgtBackground = resL$totalWgtBackground)
    percentGC <- unlist(lapply(pfms, function(x) {
        m <- TFBSTools::Matrix(x)
        100 * sum(m[c("C","G"), ]) / sum(m)
    }), use.names = FALSE)
    mid <- rep(NA, length(pfms))
    mnm <- rownames(assayL[[1]])
    if (all(grepl(":::", rownames(assayL[[1]])))) {
        mid <- sub(":::.+$", "", rownames(assayL[[1]]))
        mnm <- sub("^.+:::", "", rownames(assayL[[1]]))
        assayL <- lapply(assayL, function(x) {
            rownames(x) <- sub(":::.+$", "", rownames(x))
            x
        })
        for (i in seq_along(pfms)) {
            pfms[[i]]@ID <- mid[i]
            pfms[[i]]@name <- mnm[i]
        }
    }
    ## ... ... rowData
    rdat <- S4Vectors::DataFrame(motif.id = mid,
                                 motif.name = mnm,
                                 motif.pfm = pfms,
                                 motif.pwm = rep(NA, length(pfms)),
                                 motif.percentGC = percentGC)
    ## ... ... metadata
    mdatL <- list(regions = gr,
                  bins = b,
                  bins.binmode = attr(b, "binmode"),
                  bins.breaks = as.vector(attr(b, "breaks")),
                  bins.bin0 = getZeroBin(b),
                  param = list(method = "Homer",
                               genomedir = genomedir,
                               outdir = outdir,
                               motifFile = motifFile,
                               homerfile = homerfile,
                               regionsize = regionsize,
                               pseudocount.log2enr = pseudocount.log2enr,
                               p.adj.method = p.adjust.method,
                               Ncpu = Ncpu),
                  motif.distances = NULL)
    ## ... ... return
    SummarizedExperiment::SummarizedExperiment(
      assays = assayL, colData = cdat, rowData = rdat,
      metadata = mdatL)
}







