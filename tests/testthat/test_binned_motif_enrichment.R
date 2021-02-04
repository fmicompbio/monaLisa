context("binned motif enrichment")


test_that("get_binned_motif_enrichment() works in default mode", {

  library(GenomicRanges)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(JASPAR2018)
  library(TFBSTools)
  library(Biostrings)
  library(SummarizedExperiment)

  genome <- BSgenome.Mmusculus.UCSC.mm10

  ## We use pre-selected regions and motifs that we 
  ## ... know will be enriched in one set of regions vs the other
  ## ... Below (commented out) is a description of how we chose them:
  ##
  ## library(monaLisa)
  ## library(GenomicRanges)
  ## library(BSgenome.Mmusculus.UCSC.mm10)
  ## library(JASPAR2018)
  ## library(TFBSTools)
  ## library(Biostrings)
  ## library(SummarizedExperiment)
  ## 
  ## genome <- BSgenome.Mmusculus.UCSC.mm10
  ##  
  ## ## Find motif hits on chr10 for all TFs on 100 bp tiles from mm10
  ## 
  ## pfms <- getMatrixSet(JASPAR2018, list(matrixtype = "PFM", tax_group = "vertebrates"))
  ## pwms <- toPWM(pfms) # 579 TFs
  ##
  ## seq_info <- seqinfo(genome)
  ## seq_info <- seq_info["chr10"]
  ## gr <- tileGenome(seq_info, tilewidth = 100, cut.last.tile.in.chrom = TRUE)
  ## seqs <- getSeq(genome, gr)
  ## 
  ## hits <- findMotifHits(query = pwms, subject = seqs, min.score = 10, method = "matchPWM.concat", Ncpu = 10)
  ## 
  ## ## Select 2 TFs and the regions (bins) where they have unique hits
  ## 
  ## mat <- as.matrix(as.data.frame.matrix(table(seqnames(hits), as.character(hits$pwmname))))
  ## w <- which(mat > 1)
  ## mat_zoops <- mat
  ## mat_zoops[w] <- 1
  ## 
  ## keep <- rowSums(mat_zoops) == 1
  ## mat_zoops_unique <- mat_zoops[keep, ]
  ## dim(mat_zoops_unique)
  ## # [1] 19145   574
  ## 
  ## mat_zoops_unique[1:6, 1:6]
  ## #           Ahr::Arnt Alx1 ALX3 Alx4 Ar Arid3a
  ## # s31072         0    0    0    0  0      0
  ## # s31095         1    0    0    0  0      0
  ## # s31096         1    0    0    0  0      0
  ## # s31097         1    0    0    0  0      0
  ## # s31098         1    0    0    0  0      0
  ## # s31100         1    0    0    0  0      0
  ## 
  ## head(sort(colSums(mat_zoops_unique), decreasing = TRUE), n=14)
  ## # ZNF263                RBPJ                SPIB       Nkx2-5(var.2) SMAD2::SMAD3::SMAD4               PRDM1                 MYB 
  ## # 1055                 534                 463                 394                 346                 337                 305 
  ## # DMRT3               LIN54               Nr5a2                 VDR               SOX10                Mafb        Pou5f1::Sox2 
  ## # 286                 285                 263                 263                 257                 251                 245 
  ##
  ## ## Choose regions (bins) unique for RBPJ and SPIB respectively 
  ## RBPJ_seqnames <- rownames(mat_zoops_unique)[as.logical(mat_zoops_unique[, "RBPJ"])]
  ## SPIB_seqnames <- rownames(mat_zoops_unique)[as.logical(mat_zoops_unique[, "SPIB"])]
  ## 
  ## ## Get ranges of bins corresponding to these seqnames (seqname corresponds to index)
  ## RBPJ_ind <- as.numeric(limma::strsplit2(RBPJ_seqnames, "s")[, 2])
  ## SPIB_ind <- as.numeric(limma::strsplit2(SPIB_seqnames, "s")[, 2])
  ## 
  ## gr_RBPJ <- gr[RBPJ_ind]
  ## gr_SPIB <- gr[SPIB_ind]
  ## 
  ## gr_SPIB_final <- c(gr_SPIB[1:10], gr_RBPJ[1:2])
  ## gr_RBPJ_final <- c(gr_RBPJ[1:10], gr_SPIB[1:2])
  ##
  ## as.data.frame(gr_SPIB_final)
  ## #      seqnames   start     end width strand
  ## #   1     chr10 3133401 3133500   100      *
  ## #   2     chr10 3447301 3447400   100      *
  ## #   3     chr10 3927801 3927900   100      *
  ## #   4     chr10 4110201 4110300   100      *
  ## #   5     chr10 4179501 4179600   100      *
  ## #   6     chr10 4444901 4445000   100      *
  ## #   7     chr10 4731101 4731200   100      *
  ## #   8     chr10 5110901 5111000   100      *
  ## #   9     chr10 5114401 5114500   100      *
  ## #   10    chr10 5332501 5332600   100      *
  ## #   11    chr10 3289801 3289900   100      *
  ## #   12    chr10 3420201 3420300   100      *
  ## 
  ## as.data.frame(gr_RBPJ_final)
  ## #      seqnames   start     end width strand
  ## #   1     chr10 3289801 3289900   100      *
  ## #   2     chr10 3420201 3420300   100      *
  ## #   3     chr10 3594201 3594300   100      *
  ## #   4     chr10 3622701 3622800   100      *
  ## #   5     chr10 3877101 3877200   100      *
  ## #   6     chr10 3993201 3993300   100      *
  ## #   7     chr10 4068701 4068800   100      *
  ## #   8     chr10 4485601 4485700   100      *
  ## #   9     chr10 4641401 4641500   100      *
  ## #   10    chr10 4683701 4683800   100      *
  ## #   11    chr10 3133401 3133500   100      *
  ## #   12    chr10 3447301 3447400   100      *
  
  ## Create dataset
  gr_SPIB <- GRanges(seqnames = "chr10", 
                     ranges = IRanges(start = c(3133401, 3447301, 3927801, 4110201, 4179501, 4444901, 4731101, 5110901, 5114401, 5332501, 3289801, 3420201), 
                                      end = c(3133500, 3447400, 3927900, 4110300, 4179600, 4445000, 4731200, 5111000, 5114500, 5332600, 3289900, 3420300)), 
                     strand = "*")
  gr_RBPJ <- GRanges(seqnames = "chr10", 
                     ranges = IRanges(start = c(3289801, 3420201, 3594201, 3622701, 3877101, 3993201, 4068701, 4485601, 4641401, 4683701, 3133401, 3447301), 
                                      end = c(3289900, 3420300, 3594300, 3622800, 3877200, 3993300, 4068800, 4485700, 4641500, 4683800, 3133500, 3447400)), 
                     strand = "*")
  gr <- c(gr_RBPJ, gr_SPIB)
  names(gr) <- paste0("peak_", 1:length(gr))
  seqs <- getSeq(genome, gr)
  
  ## PWMs
  pfms <- getMatrixSet(JASPAR2018, list(matrixtype = "PFM", tax_group = "vertebrates", name = c("RBPJ", "SPIB")))
  pwms <- toPWM(pfms)
  
  ## bins
  b <- monaLisa::bin(x = 1:length(seqs),
                     binmode = "breaks",
                     breaks = c(0, (length(gr_RBPJ)), (length(seqs)+1)))
  
  ## Binned motif enrichment with monaLisa
  enr_res <- monaLisa::get_binned_motif_enrichment(seqs = seqs,
                                                   bins = b,
                                                   pwmL = pwms,
                                                   genome = genome)

  ## Tests
  # ... expected results on dataset
  expect_true(base::inherits(enr_res, "SummarizedExperiment"))
  expect_true(all(rownames(enr_res) == c("MA0081.1", "MA1116.1")))
  expect_equal(SummarizedExperiment::assay(enr_res, 1)$fg_weight_sum[1], 2)
  expect_equal(SummarizedExperiment::assay(enr_res, 1)$fg_weight_sum[2], 10)
  expect_equal(SummarizedExperiment::assay(enr_res, 1)$fg_weight_sum_total[1], 12)
  # ... missing arguments or wrong classes
  expect_error(get_binned_motif_enrichment(bins = b, pwmL = pwms, genome = genome))
  expect_error(get_binned_motif_enrichment(seqs = as.character(seqs), bins = b, pwmL = pwms, genome = genome))
  expect_error(get_binned_motif_enrichment(seqs = seqs, pwmL = pwms, genome = genome))
  expect_error(get_binned_motif_enrichment(seqs = seqs, bins = as.character(b), pwmL = pwms, genome = genome))
  expect_error(get_binned_motif_enrichment(seqs = seqs, bins = b, genome = genome))
  expect_error(get_binned_motif_enrichment(seqs = seqs, bins = b, pwmL = as.matrix(pwms[[1]]), genome = genome))
  expect_error(get_binned_motif_enrichment(seqs = seqs, bins = b, pwmL = pwms))
  expect_error(get_binned_motif_enrichment(seqs = seqs, bins = b, pwmL = pwms, genome = "mm10"))
 
})

