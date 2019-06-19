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
#' @return Absolute path to \code{homerfile}, or \code{NA} if none or several were found.
#'
#' @export
findHomer <- function(homerfile = "findMotifsGenome.pl", dirs = NULL) {
    if (is.null(dirs))
        dirs <- strsplit(x = Sys.getenv("PATH"), split = ":")[[1]]
    res <- list.files(path = dirs, pattern = paste0("^",homerfile,"$"), full.names = TRUE)
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
#' @param relScoreCutoff Currently ignored. numeric(1) in [0,1] that sets the default motif
#'     log-odds score cutof to relScoreCutoff * maximal score for each PWM
#'     (default: 0.8).
#'
#' @return \code{TRUE} if successful.
#'
#' @seealso \code{\link[TFBSTools]{getMatrixSet}} for details on the argument \code{opts}.
#'
#' @export
dumpJaspar <- function(filename, pkg = "JASPAR2018", opts = list(tax_group = "vertebrates"),
                       relScoreCutoff = 0.8) {
    stopifnot(!file.exists(filename))
    stopifnot(is.numeric(relScoreCutoff) && length(relScoreCutoff) == 1 &&
              relScoreCutoff >= 0.0 && relScoreCutoff <= 1.0)

    requireNamespace(pkg)
    requireNamespace("TFBSTools")

    # load PFMs and convert to PWMs
    mdb <- getFromNamespace(pkg, ns=pkg)
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
                           paste(TFBSTools::tags(siteList[[i]])$acc, collapse = "::"),
                           TFBSTools::tags(siteList[[i]])$type), collapse = "|")
        #the -10 is added so that the motif file has 4 columns, which is need to run compareMotifs.pl
        #for the weight matrix clustering
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
        b <- factor(b, levels=unique(b))
    stopifnot(length(b) == length(gr))
    stopifnot(is.character(outdir) && length(outdir) == 1L)
    stopifnot(is.character(motifFile) && length(motifFile) == 1L && file.exists(motifFile))
    stopifnot(is.character(homerfile) && length(homerfile) == 1L && file.exists(homerfile))

    if (file.exists(outdir))
        stop(outdir," already exists - will not overwrite existing folder")
    dir.create(outdir)

    homerFile <- file.path(outdir, "run.sh")
    fh <- file(homerFile, "w")

    message("creating foreground/background region files for HOMER")
    for(i in 1:nlevels(b)) {
        bn <- levels(b)[i]
        message("  bin ",bn)

        fgfile  <- sprintf("%s/bin_%03d_foreground.tab", outdir, i)
        bgfile  <- sprintf("%s/bin_%03d_background.tab", outdir, i)
        outputf <- sprintf("%s/bin_%03d_output", outdir, i)

        tmp.gr <- gr[b == bn]
        write.table(file=fgfile,
                    data.frame(1:length(tmp.gr), as.character(seqnames(tmp.gr)), start(tmp.gr)-1, end(tmp.gr), "+"),
                    sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

        tmp.gr <- gr[b != bn]
        write.table(file=bgfile,
                    data.frame(1:length(tmp.gr), as.character(seqnames(tmp.gr)), start(tmp.gr)-1, end(tmp.gr), "+"),
                    sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

        cat(sprintf("%s %s %s %s -bg %s -nomotif -p %d -size %s -mknown %s\n",
                    homerfile, fgfile, genomedir, outputf, bgfile, Ncpu, as.character(regionsize), motifFile), file=fh, append=TRUE)
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
    names(tabL) <- if(is.null(names(infiles))) infiles else names(infiles)
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
    fdr <- matrix(-log10(p.adjust(tmp, method="BH")), nrow=nrow(P))
    dimnames(fdr) <- dimnames(P)

    fdr[which(fdr == Inf, arr.ind = TRUE)] <- max(fdr[is.finite(fdr)])

    return(list(p=P, FDR=fdr, enr=enrTF, log2enr=log2enr))
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
#'     into bins (typically a factor).
#' @param genomedir Directory containing sequence files in Fasta format (one per chromosome).
#' @param outdir A path specifying the folder into which the output files will be written.
#' @param motifFile A file with HOMER formatted PWMs to be used in the enrichment analysis.
#' @param homerfile Path and file name of the \code{findMotifsGenome.pl} HOMER script.
#' @param regionsize The peak size to use in HOMER (\code{"given"} keeps the coordinate
#'     region, an integer value will keep only that many bases in the region center).
#' @param Ncpu Number of parallel threads that HOMER can use.
#'
#' @seealso The functions that are wrapped: \code{\link{prepareHomer}},
#'     \code{\link[base]{system}} and \code{\link{parseHomerOutput}}
#'
#' @return A list of four components (\code{p}, \code{FDR}, \code{enr} and \code{log2enr}),
#'     containing each a motif (rows) by bin (columns) matrix with raw
#'     -log10 P values, -log10 false discovery rates and motif enrichments as
#'     Pearson residuals (\code{enr}) and as log2 ratios (\code{log2enr}).
#'
#' @export
runHomer <- function(gr, b, genomedir, outdir, motifFile, homerfile = findHomer(), regionsize = "given", Ncpu=2L) {
    ## ... prepare
    message("\npreparing input files...")
    runfile <- prepareHomer(gr = gr, b = b, genomedir = genomedir, outdir = outdir,
                            motifFile = motifFile, homerfile = homerfile,
                            regionsize = regionsize, Ncpu = Ncpu)

    ## ... run
    message("\nrunning HOMER...")
    system2(command = "sh", args = runfile, env = paste0("PATH=",dirname(homerfile),":",Sys.getenv("PATH"),";"))

    ## ... parse output
    resfiles <- sprintf("%s/bin_%03d_output/knownResults.txt", outdir, seq_along(levels(b)))
    names(resfiles) <- levels(b)
    parseHomerOutput(resfiles)
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
clusterPWMs <- function(motifFile, homerdir, outfile){
  
    stopifnot(is.character(motifFile) && length(motifFile) == 1L && file.exists(motifFile))
    stopifnot(is.character(homerdir) && length(homerdir) == 1L && file.exists(homerdir))
    homerfile = findHomer("compareMotifs.pl", dirs = homerdir)
    #run
    message("running compareMotifs.pl...")
    system(sprintf("%s %s test -matrix %s", homerfile, motifFile, outfile), intern=TRUE)
    #/work/gbioinfo/Appz/Homer/Homer-4.8/bin/compareMotifs.pl test.motif test -matrix testmat
    as.matrix(read.delim(outfile, row.names=1))
}







