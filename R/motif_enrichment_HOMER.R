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
#' @description Get PWMs from the \code{JASPAR2016} or \code{JASPAR2018} package
#'     and write them into a HOMER-compatible motif file.
#'
#' @param filename Name of the output file to be created.
#' @param pkg Name of the Jaspar package to use (default: \code{JASPAR2018}).
#' @param opts a search options list used in \code{\link[TFBSTools]{getMatrixSet}}.
#' @param relScoreCutoff Currently ignored. numeric(1) in [0,1] that sets the default motif
#'     log-odds score cutof to relScoreCutoff * maximal score for each PWM
#'     (default: 0.8).
#'
#' @return \code{TRUE} if successful.
#'
#' @seealso \code{\link[TFBSTools]{getMatrixSet}} for details on the argument \code{opts}.
#'     \code{\link{homerToPFMatrixList}} to read a file with HOMER-formatted
#'     moitfs into a \code{\link[TFBSTools]{PFMatrixList}}.
#'
#' @importFrom utils getFromNamespace
#' @importFrom TFBSTools getMatrixSet Matrix name tags
#'
#' @export
dumpJaspar <- function(filename, pkg = "JASPAR2018",
                       opts = list(tax_group = "vertebrates"),
                       relScoreCutoff = 0.8) {
    .assertScalar(x = filename, type = "character")
    stopifnot(!file.exists(filename))
    .assertScalar(x = relScoreCutoff, type = "numeric", rngIncl = c(0, 1))
    if ("matrixtype" %in% names(opts) && opts[["matrixtype"]] != "PFM") {
      stop("opts[['matrixtype']] must be set to 'PFM'")
    } else {
      opts[["matrixtype"]] <- "PFM"
    }

    # load PFMs and convert to HOMER format (base probabilties)
    requireNamespace(pkg)
    mdb <- utils::getFromNamespace(pkg, ns = pkg)
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
        #scorecut <- relScoreCutoff * sum(log(apply(pwm, 2, max) / 0.25))
        scorecut <- log(2**10) # use constant cutoff for all motifs (ignore relScoreCutoff)
        pwm <- apply(pwm, 2, function(x){sprintf("%.3f", x)})
        rownames(pwm) <- tmp.rn
        wm.name <- paste(c(TFBSTools::name(siteList[[i]]),
                           paste(TFBSTools::tags(siteList[[i]])$acc, collapse = "."),
                           TFBSTools::tags(siteList[[i]])$type), collapse = "_")
        wm.name <- gsub("::", ".", wm.name)
        wm.name <- gsub("-", "", wm.name)
        wm.name <- gsub("\\s+", "", wm.name)
        if (length(grep("\\(var\\.\\d\\)", wm.name)) > 0) {
          tmp.num <- gsub("\\S+\\(var\\.(\\d)\\)\\S+", "\\1", wm.name)
          wm.name <- gsub(sprintf("\\(var\\.%s\\)", tmp.num),
                          sprintf("var%s", tmp.num), wm.name)
        }
        wm.name <- gsub("\\(PBM\\)", "", wm.name)
        wm.name.conv <- make.names(wm.name)
        if (!(wm.name == wm.name.conv)) {
          warning("Weight matrix name ", wm.name, " converted to ", wm.name.conv)
        }
        wm.name <- wm.name.conv

        #the -10 is added so that the motif file has 4 columns, which is needed to run compareMotifs.pl
        #for weight matrix clustering (4th column not used, bug in compareMotifs.pl, I think)
        cat(sprintf(">%s\t%s\t%.2f\t-10\n",
                    paste(apply(pwm, 2, function(x) { rownames(pwm)[which.max(x)] }), collapse = ""),
                    wm.name, scorecut),  file = fh, append = TRUE)
        write.table(file = fh, t(pwm), row.names = FALSE, col.names = FALSE,
                    sep = "\t", quote = FALSE, append = TRUE)
        flush(fh)
    }
    close(fh)
    message("done")

    return(TRUE)
}


#' @title Read a HOMER motif file and create a TFMatrixList.
#'
#' @description Read motifs from a file in HOMER format and create
#'     a TFMatrixList from them.
#'
#' @param filename Name of the input file with HOMER-formatted motifs.
#' @param n The number of observations (multiplied with base frequencies to
#'     create the number of observed bases at each position).
#'
#' @return A \code{\link[TFBSTools]{PFMatrixList}} with motifs from the file.
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
prepareHomer <- function(gr, b, genomedir, outdir, motifFile, homerfile = findHomer(), regionsize = "given", Ncpu=2) {
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

    if (file.exists(outdir))
        stop(outdir," already exists - will not overwrite existing folder")
    dir.create(outdir)

    homerFile <- file.path(outdir, "run.sh")
    fh <- file(homerFile, "w")

    message("creating foreground/background region files for HOMER")
    for (i in 1:nlevels(b)) {
        bn <- levels(b)[i]
        message("  bin ",bn)

        fgfile  <- sprintf("%s/bin_%03d_foreground.tab", outdir, i)
        bgfile  <- sprintf("%s/bin_%03d_background.tab", outdir, i)
        outputf <- sprintf("%s/bin_%03d_output", outdir, i)

        tmp.gr <- gr[b == bn]
        write.table(file = fgfile,
                    data.frame(1:length(tmp.gr), as.character(seqnames(tmp.gr)),
                               start(tmp.gr) - 1, end(tmp.gr), "+"),
                    sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

        tmp.gr <- gr[b != bn]
        write.table(file = bgfile,
                    data.frame(1:length(tmp.gr), as.character(seqnames(tmp.gr)),
                               start(tmp.gr) - 1, end(tmp.gr), "+"),
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
#'
#' @return A list of four components (\code{p}, \code{FDR}, \code{enr} and \code{log2enr}),
#'     containing each a motif (rows) by bin (columns) matrix with raw
#'     -log10 P values, -log10 false discovery rates and motif enrichments as
#'     Pearson residuals (\code{enr}) and as log2 ratios (\code{log2enr}).
#'
#' @export
parseHomerOutput <- function(infiles) {
    stopifnot(all(file.exists(infiles)))

    tabL <- lapply(infiles, read.delim)
    names(tabL) <- if (is.null(names(infiles))) infiles else names(infiles)
    P <- lapply(tabL, function(tab) {
        D <- tab[, c(1, 4)]
        logpVals <- -D[, 2]/log(10)
        names(logpVals) <- D[, 1]
        logpVals[order(names(logpVals))]
    })
    enrTF <- lapply(tabL, function(tab) {
        D <- tab[, c(1, 6, 7, 9)] # name, no. and % of target seqs with motif, % of bg seqs with motif
        D[, 3] <- as.numeric(sub("%$","", D[, 3])) / 100
        D[, 4] <- as.numeric(sub("%$","", D[, 4])) / 100
        obsTF <- D[, 2]
        expTF <- D[, 2] / (D[, 3] + 0.001) * (D[, 4] + 0.001)
        enr <- (obsTF - expTF) / sqrt(expTF)
        enr[ is.na(enr) ] <- 0
        names(enr) <- D[, 1]
        enr[order(names(enr))]
    })
    log2enr <- lapply(tabL, function(tab){
    	D <- tab[, c(6, 8)] # number of target seqs and bg seqs with motif
    	nTot <- as.numeric(gsub("\\S+\\.(\\d+)\\.", "\\1", colnames(D))) #total number of target and background sequences
    	D.norm <- t(min(nTot)*t(D)/nTot) # scale to smaller number (usually number of target sequences)
    	DL <- log2(D.norm + 8)
    	log2enr <- DL[, 1] - DL[, 2]
    	names(log2enr) <- tab[, 1]
    	log2enr[order(names(log2enr))]
    })

    P <- do.call(cbind, P)
    enrTF <- do.call(cbind, enrTF)
    log2enr <- do.call(cbind, log2enr)
    tmp <-  as.vector(10**(-P))
    fdr <- matrix(-log10(p.adjust(tmp, method = "BH")), nrow = nrow(P))
    dimnames(fdr) <- dimnames(P)

    fdr[which(fdr == Inf, arr.ind = TRUE)] <- max(fdr[is.finite(fdr)])

    return(list(p = P, FDR = fdr, enr = enrTF, log2enr = log2enr))
}

# internal function:  Check if all the HOMER output files already exist and if the run was successful.
# - needs the motifFile and outdir inputs to runHomer.
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
      all(sapply(motifs_outdir_list, function(x){all(motifs_motifFile %in% x)})) & (length(out_files) == nbins)

}

#' @title Prepare and run HOMER motif enrichment analysis.
#'
#' @description Run complete HOMER motif enrichment analysis, consisting of
#'     calls to \code{\link{prepareHomer}}, \code{\link[base]{system}} and
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
#' @param Ncpu Number of parallel threads that HOMER can use.
#'
#' @seealso The functions that are wrapped: \code{\link{prepareHomer}},
#'     \code{\link[base]{system}} and \code{\link{parseHomerOutput}},
#'     \code{\link{bin}} for binning of regions
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} \code{x}
#'   with
#'   \describe{
#'     \item{assays(x)}{containing the four components \code{p}, \code{FDR},
#'         \code{enr} and \code{log2enr}), each a motif (rows) by bin (columns)
#'         matrix with raw -log10 P values, -log10 false discovery rates and motif
#'         enrichments as Pearson residuals (\code{enr}) and as log2 ratios
#'         (\code{log2enr}).}
#'     \item{rowData(x)}{containing information about the motifs.}
#'     \item{colData(x)}{containing information about the bins.}
#'     \item{metaData(x)}{containing meta data on the object (e.g. parameter
#'         values and distances between pairs of motifs).}
#'   }
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#'
#' @export
runHomer <- function(gr, b, genomedir, outdir, motifFile,
                     homerfile = findHomer(), regionsize = "given", Ncpu = 2L) {
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
    .assertScalar(x = Ncpu, type = "numeric", rngIncl = c(1, Inf))

    ## ... check if the HOMER output is already there for all bins and if it ran completely:
    ## ... ... If yes, go to the 'parse output step', otherwise run homer and check again
    if (.checkHomerRun(motifFile = motifFile, outdir = outdir, nbins = nlevels(b))) {

      message("\nHOMER output files already exist, using existing files...")

    } else {

      ## ... case: all/some files exist and/or HOMER didn't run correctly: warn the User to remove all existing files and rerun
      if (!.checkHomerRun(motifFile = motifFile, outdir = outdir, nbins = nlevels(b)) &&
          any(file.exists(dir(path = outdir, pattern = "knownResults.txt",
                              full.names = TRUE, recursive = TRUE, ignore.case = FALSE)))) {

        stop("\nThere are existing 'knownResults.txt' file(s) in outdir. ",
             "There may be missing 'knownResults.txt' files for some bins ",
             "and/or the existing files are incomplete (cases where the HOMER ",
             "run failed). Please delete these files and rerun 'runHomer'.")

      }

      ## ... prepare
      message("\npreparing input files...")
      runfile <- prepareHomer(gr = gr, b = b, genomedir = genomedir, outdir = outdir,
                              motifFile = motifFile, homerfile = homerfile,
                              regionsize = regionsize, Ncpu = Ncpu)

      ## ... run
      message("\nrunning HOMER...")
      system2(command = "sh", args = runfile, env = paste0("PATH=",dirname(homerfile),":",Sys.getenv("PATH"),";"))

      ## ... check HOMER ran correctly
      if (!.checkHomerRun(motifFile = motifFile, outdir = outdir, nbins = length(levels(b)))) {
          stop("HOMER output wasn't complete. Try running again.")
      }

    }

    ## ... parse output
    resfiles <- sprintf("%s/bin_%03d_output/knownResults.txt", outdir, seq_along(levels(b)))
    names(resfiles) <- levels(b)
    resL <- parseHomerOutput(resfiles)

    ## ... create SummarizedExperiment
    cdat <- S4Vectors::DataFrame(bin.names = levels(b),
                                 bin.lower = attr(b, "breaks")[-(nlevels(b) + 1)],
                                 bin.upper = attr(b, "breaks")[-1],
                                 bin.nochange = seq.int(nlevels(b)) %in% attr(b, "bin0"))
    cdat <- cdat[match(colnames(resL[[1]]), levels(b)), ]
    pfms <- homerToPFMatrixList(motifFile)
    pfms <- pfms[match(rownames(resL[[1]]), TFBSTools::name(pfms))]
    percentGC <- unlist(lapply(pfms, function(x) {
        m <- TFBSTools::Matrix(x)
        100 * sum(m[c("C","G"), ]) / sum(m)
    }), use.names = FALSE)
    rdat <- S4Vectors::DataFrame(motif.name = rownames(resL[[1]]),
                                 motif.pfm = pfms,
                                 motif.percentGC = percentGC)
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = resL, colData = cdat, rowData = rdat,
      metadata = list(regions = gr,
                      bins = b,
                      bins.binmode = attr(b, "binmode"),
                      bins.breaks = as.vector(attr(b, "breaks")),
                      bins.bin0 = attr(b, "bin0"),
                      param.genomedir = genomedir,
                      param.outdir = outdir,
                      param.motifFile = motifFile,
                      param.homerfile = homerfile,
                      param.regionsize = regionsize,
                      param.Ncpu = Ncpu,
                      motif.distances = NULL)
    )
    return(se)
}







