context("QualityModString")
test_that("QualityModString:",{
  dnaTestSeq <- paste(alphabet(ModDNAString()), collapse = "")
  seq <- ModDNAString(dnaTestSeq)
  qseq <- PhredQuality(paste0(rep("!",
                                  length(seq)),
                              collapse = ""))
  qset <- QualityScaledModDNAStringSet(seq,
                                       qseq)
  expect_equal(as.character(qset),as.character(seq))
  qset <- QualityScaledModDNAStringSet(list(seq,seq),
                                       c(qseq,qseq))
  expect_equal(as.character(qset),as.character(ModDNAStringSet(list(seq,seq))))
  #
  rnaTestSeq <- paste(alphabet(ModRNAString()), collapse = "")
  seq <- ModRNAString(rnaTestSeq)
  qseq <- PhredQuality(paste0(rep("!",
                                  length(seq)),
                              collapse = ""))
  qset <- QualityScaledModRNAStringSet(seq,
                                       qseq)
  expect_equal(as.character(qset),as.character(seq))
  qset <- QualityScaledModRNAStringSet(list(seq,seq),
                                       c(qseq,qseq))
  expect_equal(as.character(qset),as.character(ModRNAStringSet(list(seq,seq))))
})
