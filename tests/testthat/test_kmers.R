context("k-mers")

test_that(".cons2matrix works as expected", {
    expect_error(.cons2matrix(x = c("error","error")))
    expect_error(.cons2matrix(x = 1L))
    expect_error(.cons2matrix("ACGT", n = "error"))
    
    expect_is(res1 <- .cons2matrix(x = "ACGT", n = 1L), "matrix")
    expect_is(res2 <- .cons2matrix(x = "ACGT", n = 2L), "matrix")
    expect_identical(dim(res1), c(4L, 4L))
    expect_identical(rownames(res1), c("A","C","G","T"))
    expect_equal(res1, diag(4), check.attributes = FALSE)
    expect_identical(res1 * 2L, res2)
})

test_that("getKmerFreq works as expected", {
    library(Biostrings)

    ## truly random sequences...
    set.seed(1)
    seqs <- sapply(1:100, function(i) paste(sample(x = c("A","C","G","T"),
                                                   size = 500L, replace = TRUE),
                                            collapse = ""))
    r <- sample(500L - 5L, 100L)
    ## ... with a planted 6-mer
    substr(seqs, start = r, stop = r + 5L) <- "AACGTT"
    seqsDSS <- Biostrings::DNAStringSet(seqs)

    expect_error(getKmerFreq(1L))
    expect_error(getKmerFreq(seqs, kmerLen = "error"))
    expect_error(getKmerFreq(seqs, kmerLen = 1:2))
    expect_error(getKmerFreq(seqs, kmerLen = 2.5))
    expect_error(getKmerFreq(seqsDSS, MMorder = "error"))
    expect_error(getKmerFreq(seqsDSS, MMorder = 1:2))
    expect_error(getKmerFreq(seqsDSS, MMorder = 2.5))
    expect_error(getKmerFreq(seqsDSS, MMorder = -1))
    expect_error(getKmerFreq(seqsDSS, MMorder = 4))
    expect_error(getKmerFreq(seqs, 3, includeRevComp = "error"))

    ## zoops = FALSE
    expect_is(res1 <- getKmerFreq(seqs, kmerLen = 4, zoops = FALSE),    "list")
    expect_is(res2 <- getKmerFreq(seqsDSS, kmerLen = 4, zoops = FALSE), "list")
    expect_identical(res1, res2)
    expect_equal(sum(res1$freq.obs), sum(res1$freq.exp), tolerance = 0.001)
    expect_identical(names(res1$log2enr)[which.max(res1$log2enr)], "ACGT")
    expect_true(all(names(res1$padj)[res1$padj < 0.001] %in% c("AACG", "ACGT", "CGTT")))

    ## zoops = TRUE
    set.seed(2)
    seqs <- sapply(1:100, function(i) paste(sample(x = c("A","C","G","T"),
                                                   size = 500L, replace = TRUE),
                                            collapse = ""))
    seqs <- c(DNAStringSet(DNAString(paste(rep("CA", 250), collapse = ""))), seqs[-1])
    expect_is(res3 <- getKmerFreq(seqs, kmerLen = 4, zoops = TRUE),  "list")
    expect_is(res4 <- getKmerFreq(seqs, kmerLen = 4, zoops = FALSE, includeRevComp = FALSE), "list")
    expect_equal(sum(res3$padj < 0.001), 0L)
    expect_equal(sum(res4$padj < 0.001), 2L)
    expect_equal(names(res4$padj[res4$padj < 0.001]), c("ACAC", "CACA"))

    ## strata
    expect_is(res5 <- getKmerFreq(seqsDSS, kmerLen = 4, zoops = FALSE, strata = 3), "list")
    expect_equal(sum(res1$freq.obs), sum(res5$freq.obs), tolerance = 0.001)
    expect_equal(sum(res1$freq.exp), sum(res5$freq.exp), tolerance = 0.001)
    expect_true(all(res5$strata %in% c(1L, 2L, 3L)))
    expect_length(res5$freq.strata, 3L)
    expect_equal(sum(res5$CpGoe), 204.030341557169, tolerance = 1e-6)
    expect_identical(names(res5$log2enr)[which.max(res5$log2enr)], "ACGT")
    expect_true(all(names(res5$padj)[res5$padj < 0.001] %in% c("AACG", "ACGT", "CGTT")))
})

test_that("countKmerPairs and countKmerPairsSelected work as expected", {

    ### input sequences of variable length with and without N bases
    seqs <- Biostrings::DNAStringSet(c("AAAAAAN","ATATAT","ACGTAC","N"))
    kmerschar <- c("AA","AC","AT","CG","GT","TA")
    kmers <- Biostrings::DNAStringSet(kmerschar)

    # ... countKmerPairs
    expect_error(countKmerPairs(x = "error"))
    expect_error(countKmerPairs(x = seqs, k = 0))
    expect_error(expect_warning(countKmerPairs(x = seqs, k = 20, n = 0)))
    expect_error(countKmerPairs(x = seqs, k = 2, n = 0))

    expect_is(res1 <- countKmerPairs(x = seqs, k = 2, n = 1, zoops = FALSE), "dgCMatrix")
    expect_is(res2 <- countKmerPairs(x = seqs, k = 2, n = 1, zoops = TRUE),  "dgCMatrix")
    expect_identical(res1 > 0, res2 > 0)
    expect_equal(unname(which(as.matrix(res1) > 0, arr.ind = TRUE)),
                 cbind(c(1, 13, 13, 2, 7,  4,  12),
                       c(1,  2,  4, 7, 12, 13, 13)))
    expect_identical(sum(res1), 12)
    expect_identical(sum(res2), 7)
    
    # ... countKmerPairsSelected
    expect_error(countKmerPairsSelected(x = "error", kmers = "error"))
    expect_error(countKmerPairsSelected(x = seqs, kmers = "error"))
    expect_error(countKmerPairsSelected(x = seqs, kmers = Biostrings::DNAStringSet(c("A","AA"))))
    expect_error(countKmerPairsSelected(x = seqs, kmers = Biostrings::DNAStringSet(c("A","A"))))
    expect_error(expect_warning(countKmerPairsSelected(x = seqs, kmers = Biostrings::DNAStringSet(rep("AA",60001)), n = 0)))
    expect_error(countKmerPairsSelected(x = seqs, kmers = kmers, n = 0))

    expect_is(res1b <- countKmerPairsSelected(x = seqs, kmers = kmers, n = 1, zoops = FALSE), "dgCMatrix")
    expect_is(res2b <- countKmerPairsSelected(x = seqs, kmers = kmers, n = 1, zoops = TRUE),  "dgCMatrix")
    expect_identical(res1b > 0, res2b > 0)
    expect_identical(dimnames(res1b), dimnames(res2b))
    expect_identical(dimnames(res1b), list(kmerschar, kmerschar))
    expect_identical(res1[kmerschar, kmerschar], res1b)
    expect_identical(res2[kmerschar, kmerschar], res2b)
    expect_true(sum(res1b) < sum(width(seqs) - 1))

    ### alternative input sequences
    seqs2 <- Biostrings::DNAStringSet(c("AAANAANAANAAN"))
    
    # ... countKmerPairs
    expect_is(res3 <- countKmerPairs(x = seqs2, k = 2, n = 1, zoops = TRUE),  "dgCMatrix")
    expect_is(res4 <- countKmerPairs(x = seqs2, k = 2, n = 1, zoops = FALSE), "dgCMatrix")
    expect_identical(res3[1,1], 1)
    expect_true(all(as.vector(res3[-1]) == 0))
    expect_identical(dim(res3), c(16L, 16L))
    expect_identical(res3, res4)
    
    # ... countKmerPairsSelected
    expect_is(res3b <- countKmerPairsSelected(x = seqs2, kmers = kmers, n = 1, zoops = TRUE),  "dgCMatrix")
    expect_is(res4b <- countKmerPairsSelected(x = seqs2, kmers = kmers, n = 1, zoops = FALSE), "dgCMatrix")
    expect_identical(res3b[1,1], 1)
    expect_true(all(as.vector(res3b[-1]) == 0))
    expect_identical(dim(res3b), c(length(kmers), length(kmers)))
    expect_identical(res3, res4)
})

test_that("clusterKmers works as expected", {
    library(Biostrings)

    ## truly random sequences...
    set.seed(1)
    seqs <- sapply(1:100, function(i) paste(sample(x = c("A","C","G","T"),
                                                   size = 500L, replace = TRUE),
                                            collapse = ""))
    dseqs <- Biostrings::DNAStringSet(seqs)
    r <- sample(500L - 5L, 100L)
    ## ... with a planted 6-mer
    substr(seqs, start = r, stop = r + 5L) <- "AACGTT"
    x1 <- getKmerFreq(seqs, kmerLen = 4, zoops = FALSE, includeRevComp = FALSE)
    x2 <- names(x1$padj[order(x1$padj)[1:10]])
    expect_error(clusterKmers(x2, method = "cooccurrence"))
    expect_message(res1 <- clusterKmers(x1, method = "similarity"))
    expect_message(res2 <- clusterKmers(x2, method = "similarity"))
    expect_message(res3 <- clusterKmers(x2, method = "similarity",
                                        includeRevComp = TRUE))
    expect_message(res4 <- clusterKmers(x2, method = "cooccurrence",
                                        seqs = dseqs))
    expect_message(res5 <- clusterKmers(x2, method = "cooccurrence",
                                        seqs = dseqs, includeRevComp = TRUE))
    expect_message(expect_warning(res6 <- clusterKmers(x1, method = "cooccurrence",
                                                       seqs = dseqs, n = 11)))

    expect_type(res1, "double")
    expect_length(res1, 10L)
    expect_equal(res1, res2, check.attributes = FALSE)
    expect_identical(as.vector(res1), c(1, 1, 1, 2, 3, 1, 3, 1, 4, 5))
    expect_identical(as.vector(res3), c(1, 1, 1, 2, 3, 1, 3, 1, 3, 4))
    expect_identical(as.vector(res4), c(5, 1, 5, 2, 3, 5, 4, 5, 6, 7))
    expect_identical(as.vector(res5), c(5, 1, 5, 2, 3, 5, 4, 5, 6, 7, 8, 10, 9, 10, 11, 12, 13))
    expect_identical(as.vector(res6), c(2, 3, 2, 1, 3, 2, 3, 2, 3, 2))
})

test_that(".calcKmerEnrichment works", {
    fseq <- c(f1  = "AGGGGGGGGG", f3  = "AAAGGGGGGG", f4  = "AAAAGGGGGG",
              f6  = "AAAAAAGGGG", f8  = "AAAAAAAAGG")
    bseq <- c(b1  = "AGGGGGGGGG", b2  = "AAGGGGGGGG", b3  = "AAAGGGGGGG",
              b4  = "AAAAGGGGGG", b5  = "AAAAAGGGGG")
    df <- DataFrame(seqs = DNAStringSet(c(fseq, bseq)),
                    isForeground = rep(c(TRUE, FALSE), c(length(fseq), length(bseq))),
                    GCfrac = NA_real_,
                    GCbin = NA_integer_,
                    GCwgt = NA_real_,
                    seqWgt = NA_real_)
    attr(df, "err") <- 0
    df <- .calculateGCweight(df)
    df <- .iterativeNormForKmers(df)
    k <- 2L

    expect_error(.calcKmerEnrichment(k = "error", df = df),
                 "'k' must be of type 'numeric'")
    expect_error(.calcKmerEnrichment(k = k, df = "error"),
                 "'df' should be a DataFrame")
    expect_error(.calcKmerEnrichment(k = k, df = df,
                                     test = "error"),
                 "should be one of")
    expect_error(.calcKmerEnrichment(k = k, df = df,
                                     test = "fisher", verbose = "error"),
                 "'verbose' must be of type 'logical'")
    
    expect_message(res1 <- .calcKmerEnrichment(k = k, df = df,
                                               test = "binomial", verbose = TRUE))
    expect_message(res3 <- .calcKmerEnrichment(k = k, df = df,
                                               test = "fisher", verbose = TRUE))
    
    expect_is(res1, "data.frame")
    expect_is(res3, "data.frame")
    
    expect_equal(dim(res1), c(4^k, 6))
    expect_identical(res1$sumForegroundWgtWithHits, res3$sumForegroundWgtWithHits)
    expect_identical(res1$sumBackgroundWgtWithHits, res3$sumBackgroundWgtWithHits)
    expect_equal(res1$logP, c(-0.163497930965412, 0, -0.843717559537887, 0, 0,
                              0, 0, 0, 0, 0, -0.843717559537887, 0, 0, 0, 0, 0))
    expect_equal(res3$logP, c(-0.154150679827258, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0))
})


test_that("calcBinnedKmerEnr works as expected", {
    library(SummarizedExperiment)
    library(BiocParallel)
    pparams <- MulticoreParam(2L)

    set.seed(1)
    k <- 3L
    seqstr <- unlist(lapply(1:200,
                            function(i) paste(sample(c("A","C","G","T"),
                                                     50, replace = TRUE),
                                              collapse = "")))
    b <- bin(rep(1:2, each = 100), binmode = "equalN", nElements = 100)
    m <- structure(c("AAC", "CCA"), names = levels(b))
    mrc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(m)))
    m2 <- rbind(m, mrc)
    rownames(m2) <- NULL
    for (b1 in levels(b))
        substr(seqstr[b == b1], start = 10, stop = 12) <- m[b1]
    seqs <- DNAStringSet(seqstr)
    gnm <- DNAStringSet(unlist(lapply(1:10,
                                      function(i) paste(sample(c("A","C","G","T"),
                                                               1000 - i * 10,
                                                               replace = TRUE),
                                                        collapse = ""))))
    names(gnm) <- paste0("g", seq_along(gnm))
    

    expect_error(calcBinnedKmerEnr("error"))
    expect_error(calcBinnedKmerEnr(seqs, bins = NULL))
    expect_error(calcBinnedKmerEnr(seqs, bins = b[1:20]))
    expect_error(calcBinnedKmerEnr(seqs, bins = as.numeric(b)))
    expect_error(calcBinnedKmerEnr(seqs, b, kmerLen = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, background = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, background = "model", MMorder = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, background = "zeroBin"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, background = "genome"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, test = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, includeRevComp = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, maxFracN = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, maxKmerSize = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, GCbreaks = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, pseudocount.kmers = -1))
    expect_error(calcBinnedKmerEnr(seqs, b, k, pseudocount.log2enr = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, pseudofreq.pearsonResid = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, zoops = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, p.adjust.method = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, background = "genome",
                                   genome = gnm,
                                   genome.regions = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, background = "genome",
                                   genome = gnm,
                                   genome.regions = GRanges("error",
                                                            IRanges(1, 10))))
    expect_error(calcBinnedKmerEnr(seqs, b, k, background = "genome",
                                   genome = gnm, genome.oversample = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, background = "genome",
                                   genome = gnm, genome.seed = "error"))
    expect_error(calcBinnedKmerEnr(seqs, b, k, BPPARAM = "error"))
    expect_error(calcBinnedKmerEnr(DNAStringSet(rep("NNNNNNNNNN", 10)), background = "model"))

    expect_message(res1 <- calcBinnedKmerEnr(seqs, b, k, includeRevComp = FALSE, verbose = TRUE))
    RNGkind("L'Ecuyer-CMRG")
    set.seed(42L)
    res2 <- calcBinnedKmerEnr(seqs, b, k, background = "genome", genome = gnm,
                              includeRevComp = FALSE, verbose = FALSE, BPPARAM = pparams)
    RNGkind("default")
    res3 <- calcBinnedKmerEnr(seqs, b, k, background = "model",
                              BPPARAM = pparams)
    res4 <- calcBinnedKmerEnr(seqs, b, k, background = "model",
                              test = "binomial")

    expect_is(res1, "SummarizedExperiment")
    expect_is(res2, "SummarizedExperiment")
    expect_is(res3, "SummarizedExperiment")
    expect_is(res4, "SummarizedExperiment")
    expect_identical(assayNames(res1),
                     c("negLog10P", "negLog10Padj", "pearsonResid", "log2enr",
                       "sumForegroundWgtWithHits", "sumBackgroundWgtWithHits"))
    expect_identical(apply(assay(res1, "negLog10Padj"), 2,
                           function(x) names(x)[x > 3]), m)
    expect_identical(apply(assay(res2, "negLog10Padj"), 2,
                           function(x) names(x)[x > 3]), m)
    expect_identical(apply(assay(res3, "negLog10Padj"), 2,
                           function(x) names(x)[x > 0.5]), m2)
    expect_identical(apply(assay(res4, "negLog10Padj"), 2,
                           function(x) names(x)[x > 1]), m2)
    expect_equal(colSums(assay(res1, "negLog10P")),
                 c(`[1,1.5]` = 33.1936891496806, `(1.5,2]` = 31.5395993919718))
    if(.Platform$OS.type == "unix") {
        expect_equal(colSums(assay(res2, "negLog10P")),
                     c(`[1,1.5]` = 37.1270049027622, `(1.5,2]` = 40.6593382003256))
    }
    expect_equal(colSums(assay(res3, "negLog10P")),
                 c(`[1,1.5]` = 27.2409566382244, `(1.5,2]` = 23.5971145106083))
    expect_equal(colSums(assay(res4, "negLog10P")),
                 c(`[1,1.5]` = 37.2309483817477, `(1.5,2]` = 30.1142668663514))
    expect_identical(nrow(res1), as.integer(4^k))
    expect_identical(ncol(res1), nlevels(b))
    expect_identical(assay(res1, "sumForegroundWgtWithHits"),
                     assay(res2, "sumForegroundWgtWithHits"))
    expect_identical(assays(res3)[-c(1, 2)], assays(res4)[-c(1, 2)])
})

test_that("convertKmersToMotifs works as expected", {
    library(SummarizedExperiment)
    library(TFBSTools)

    ## truly random sequences...
    nseqs <- 1000
    nbins <- 5
    seqlen <- 100
    set.seed(1)
    seqs <- sapply(seq.int(nseqs), function(i) paste(sample(x = c("A","C","G","T"),
                                                            size = seqlen, replace = TRUE),
                                                     collapse = ""))
    b <- bin(rep(seq.int(nbins), each = round(nseqs / nbins)), binmode = "equalN", nElements = round(nseqs / nbins))
    ## ... with a planted 6-mer
    for (ib in 1:3) {
        i <- which(as.numeric(b) == ib)
        i <- i[seq.int(round(length(i) / ib))]
        r <- sample(seqlen - 5L, length(i), replace = TRUE)
        substr(seqs[i], start = r, stop = r + 5L) <- "AACGTT"
    }
    seqs <- Biostrings::DNAStringSet(seqs)
    res1 <- calcBinnedKmerEnr(seqs, b, kmerLen = 4, background = "model", includeRevComp = FALSE, verbose = TRUE)
    #o <- order(assay(res1, "log2enr")[, 1], decreasing = TRUE)[1:10]
    #res2 <- plotMotifHeatmaps(res1[o, ], cluster = TRUE, show_dendrogram = TRUE)

    ## motifs...
    pfms <- do.call(PFMatrixList, list(PFMatrix(ID = "m1", name = "m1",
                                                profileMatrix = rbind(A = c(85, 85,  5,  5,  5,  5),
                                                                      C = c( 5,  5, 85,  5,  5,  5),
                                                                      G = c( 5,  5,  5, 85,  5,  5),
                                                                      T = c( 5,  5,  5,  5, 85, 85))),
                                       PFMatrix(ID = "m2", name = "m2",
                                                profileMatrix = rbind(A = c( 5,  5,  5,  5),
                                                                      C = c(85,  5,  5, 85),
                                                                      G = c( 5, 85, 85,  5),
                                                                      T = c( 5,  5,  5,  5)))))
    tf <- tempfile(fileext = ".motif")
    monaLisa:::.dumpPWMsToHomer2File(pwmL = toPWM(pfms), fname = tf)
    pfms <- homerToPFMatrixList(tf) # read back in order to have identical @listData@tags
    a1 <- convertKmersToMotifs(res1, tf, verbose = TRUE)
    a2 <- convertKmersToMotifs(res1, pfms)
    expect_is(a1, "SummarizedExperiment")
    expect_is(a2, "SummarizedExperiment")
    expect_identical(a1, a2)
    expect_equal(rowSums(assay(a1, "pearsonResid")), c(`m1:::m1` = 52.7476459257492, `m2:::m2` = -0.101202243464047))
})

test_that("extractOverlappingKmerFrequencies works as expected", {
    seqs <- Biostrings::DNAStringSet(c(s1 = "AAAAACCGTTAAAAAAAAAAAAAAAAACCGTTAAAAAAAAAAAAAA",
                                       s2 = "AAAAAAAAAAAAAAAAACCGTTAAAAAAAAAAATAACGGAAAAAAA",
                                       s3 = "ATAACGGAAAACCGTTCCGTTAAAAAAAAAAAAAAAAAAAAAAAAA"))
    seqsrc <- Biostrings::reverseComplement(seqs)
    names(seqsrc) <- paste0("r", seq_along(seqsrc))
    kmers <- c("CCGT", "CGTT", "GTTA")
    kmersrc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(kmers)))
    
    expect_error(extractOverlappingKmerFrequencies(seqs = "error", x = kmers))
    expect_error(extractOverlappingKmerFrequencies(seqs, x = "error"))
    expect_error(extractOverlappingKmerFrequencies(seqs, kmers, BPPARAM = "error"))

    res1 <- extractOverlappingKmerFrequencies(seqs, kmers)
    res2 <- extractOverlappingKmerFrequencies(seqs, kmers, BPPARAM = BiocParallel::MulticoreParam(2L))
    res3 <- extractOverlappingKmerFrequencies(c(seqs, seqsrc), kmers, includeRevComp = FALSE)
    res4 <- extractOverlappingKmerFrequencies(seqs, c(kmers, kmersrc), includeRevComp = FALSE)
    
    expect_is(res1, "integer")
    expect_identical(res1, c(CCGTTA = 5L, CCGTTCCGTTA = 1L))
    expect_identical(res1, res2)
    expect_identical(res1, res3)
    expect_identical(res4, c(CCGTTA = 3L, TAACGG = 2L, CCGTTCCGTTA = 1L))
})

test_that("buildDirGraphFromKmers works as expected", {
    set.seed(1)
    x <- sapply(1:100, function(i) paste(sample(x = c("A","C","G","T"),
                                                size = 50L, replace = TRUE),
                                         collapse = ""))
    seqs <- DNAStringSet(x = x)
    enrkmers <- c("CATAA", "GTCGA")
    g <- buildDirGraphFromKmers(seqs, enrkmers)
    
    expect_is(g, "igraph")
    expect_true(all(enrkmers %in% vertex_attr(g, "name")))
    
    olk <- extractOverlappingKmerFrequencies(seqs, enrkmers)
    onf <- oligonucleotideFrequency(DNAStringSet(names(olk)), width = 5)
    onf <- onf[, colSums(onf) > 0]
    
    expect_true(all(sort(vertex_attr(g, "name")) == sort(colnames(onf))))
})

