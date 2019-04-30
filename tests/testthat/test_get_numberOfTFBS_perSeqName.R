test_that("get_numberOfTFBS_perSeqName() works properly", {

  # create data set
  pwms <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018,
                                  list(matrixtype = "PWM", tax_group = "vertebrates",
                                       ID = c("MA0095.2", "MA0701.1", "MA0741.1")))
  sbj <- GenomicRanges::GRanges(seqnames = rep("chr1", 2), strand = "*",
                                ranges = IRanges::IRanges(start = c(9e7, 1e8), width = 5000))
  tf_gr <- findMotifHits(query = pwms, subject = sbj, min.score = 8.0, method = "matchPWM",
                         genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10, Ncpu = 1L)

  # apply function
  m <- get_numberOfTFBS_perSeqName(TFBS_gr = tf_gr, subject_gr = sbj, PWMs = pwms, Ncpu = 1L)

  # tests
  expect_true(inherits(m, "matrix"))
  expect_true(all(dim(m) == c(2,3)))


})
