context("binned motif enrichment")


test_that("get_binned_motif_enrichment() works in default mode", {

  library(GenomicRanges)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(JASPAR2018)
  library(TFBSTools)
  library(Biostrings)
  library(SummarizedExperiment)

  genome <- BSgenome.Mmusculus.UCSC.mm10

  ## Create dataset
  # ... GRanges and seqs of foreGround and Background
  gr_fg <- GRanges(seqnames =c("chr12"),
                   ranges = IRanges(start = seq(from = 85473000, by = 500, length.out = 10), width = 200),
                   strand = "*")
  seqlevels(gr_fg) <- c("chr12", "chr18")
  gr_bg <- GRanges(seqnames =c("chr18"),
                   ranges = IRanges(start = seq(from = 34850000, by = 500, length.out = 10), width = 200),
                   strand = "*")
  seqlevels(gr_bg) <- c("chr12", "chr18")
  gr <- c(gr_fg, gr_bg)
  names(gr) <- paste0("peak_", 1:length(gr))
  seqs <- getSeq(genome, gr)
  # ... bin into foreGround and backGround
  b <- bin(x = 1:length(seqs),
           binmode = "breaks",
           breaks = c(0, (length(gr_fg)), (length(seqs)+1)))
  # ... get PWMs list
  pfms <- getMatrixSet(JASPAR2018, list(matrixtype = "PFM", tax_group = "vertebrates", ID = c("MA0528.1")))
  pwms <- toPWM(pfms)

  ## Binned motif enrichment
  enr_res <- get_binned_motif_enrichment(seqs = seqs,
                                         bins = b,
                                         pwmL = pwms,
                                         genome = genome)
  
  ## Tests
  # expect_type(enr_res, "SummarizedExperiment")
  expect_true(base::inherits(enr_res, "SummarizedExperiment"))
  expect_equal(SummarizedExperiment::assay(enr_res, 1)$motif_name[1], "MA0528.1")
  expect_equal(SummarizedExperiment::assay(enr_res, 1)$fg_weight_sum[1], 2)
  expect_equal(SummarizedExperiment::assay(enr_res, 1)$fg_weight_sum_total[1], 5)
  
  # ... missing BSgenome or wrong class
  expect_error(get_binned_motif_enrichment(seqs = seqs, bins = b, pwmL = pwms))
  expect_error(get_binned_motif_enrichment(seqs = seqs, bins = b, pwmL = pwms, genome = "mm10"))
 
  # seqs <- getSeq(genome)
  # hits <- findMotifHits(query = pwms, subject = seqs[[12]], min.score = 10, method = "matchPWM.concat", Ncpu = 10, genome=genome)
  # 
  # 
  # # we take the top 6 hits for MA0528.1 and MA0476.1 on chr12, and flank by 300 bp (they are non-overlapping)
  # gr_fg <- GRanges(seqnames =c("chr12"),
  #                  ranges = IRanges(start = c(8088280, 26542048, 28904640, 49914565, 51868008, 52411520), 
  #                                   end = c(8088300, 26542068, 28904660, 49914585, 51868028, 52411540)),
  #                  strand = "*")
  # gr_fg <- gr_fg + 500
  # 
  # gr_bg <- GRanges(seqnames =c("chr12"),
  #                  ranges = IRanges(start = c(6243216, 6964571, 7569971, 7996102, 8142629, 8400584), 
  #                                   end = c(6243226, 6964581, 7569981, 7996112, 8142639, 8400594)),
  #                  strand = "*")
  # gr_bg <- gr_bg + 500
  # 
  # gr <- c(gr_fg, gr_bg)
  # names(gr) <- paste0("peak_", 1:length(gr))
  # seqs <- getSeq(genome, gr)
  # 
  # # ... bin into foreGround and backGround
  # b <- bin(x = 1:length(seqs),
  #          binmode = "breaks",
  #          breaks = c(0, (length(seqs)/2), (length(seqs)+1)))
  # 
  # # ... get PWMs list
  # pfms <- getMatrixSet(JASPAR2018, list(matrixtype = "PFM", tax_group = "vertebrates", ID = c("MA0528.1", "MA0476.1")))
  # pwms <- toPWM(pfms)
  # 
  # ## Binned motif enrichment
  # enr_res <- get_binned_motif_enrichment(seqs = seqs,
  #                                        bins = b,
  #                                        pwmL = pwms,
  #                                        genome = genome)
  # 
  # 
  # ## take fg = fos gene regions
  # start <- seq(from = 85473000, to = 85477273, by = 200)
  # gr_fg <- GRanges(seqnames ="chr12",
  #                  ranges = IRanges(start = start[1:(length(start)-1)], 
  #                                   end = start[2:length(start)]-1),
  #                  strand = "*")
  # ## take bg = egr1 gene regions
  # start <- seq(from = 34859000, to = 34864984, by = 200)
  # gr_bg <- GRanges(seqnames ="chr18",
  #                  ranges = IRanges(start = start[1:(length(start)-1)], 
  #                                   end = start[2:length(start)]-1),
  #                  strand = "*")
  # 
  # gr <- c(gr_fg, gr_bg)
  # names(gr) <- paste0("peak_", 1:length(gr))
  # seqs <- getSeq(genome, gr)
  # 
  # ## look for enrichment of fos motif
  # # ... bin into foreGround and backGround
  # b <- bin(x = 1:length(seqs),
  #          binmode = "breaks",
  #          breaks = c(0, (length(gr_fg)), (length(seqs)+1)))
  # 
  # # ... get PWMs list
  # pfms <- getMatrixSet(JASPAR2018, list(matrixtype = "PFM", tax_group = "vertebrates", ID = c("MA0476.1")))
  # pwms <- toPWM(pfms)
  # # ... binned motif enrichment
  # enr_res <- get_binned_motif_enrichment(seqs = seqs,
  #                                        bins = b,
  #                                        pwmL = pwms,
  #                                        genome = genome)
  # 
  # findMotifHits(query = pwms, subject = seqs, min.score = 10, method = "matchPWM.concat", Ncpu = 10, genome = genome)

})




