test_that("get_numberOfTFBS_perSeqName() works properly", {
  
  # create data set
  pwms <- getMatrixSet(JASPAR2018, list(matrixtype="PWM", tax_group="vertebrates"))[c(120,250,300)]
  sbj <- GRanges(seqnames = rep("chr1", 2), strand = "*", ranges = IRanges(start = c(90000000, 100000000), width = 5000))                                     
  tf_gr <- findMotifHits(query = pwms, subject = sbj, min.score = 8.0, method = "matchPWM", genome = genome, Ncpu = 1L)
  
  # apply function
  m <- get_numberOfTFBS_perSeqName(TFBS_gr = tf_gr, subject_gr = sbj, PWMs = pwms, nCpu = 1L) 
  
  # tests
  expect_true(inherits(m, "matrix"))
  expect_true(all(dim(m)==c(2,3)))
 
  
})
