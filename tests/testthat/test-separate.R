context("ModString separate/combine")
test_that("ModString separate/combine:",{
  dnaTestSeq <- paste(alphabet(ModDNAString()), collapse = "")
  seq <- ModDNAString(dnaTestSeq)
  gr <- separate(seq)
  expect_equal(length(gr),47)
  seq2 <- combineIntoModstrings(as(seq,"DNAString"),gr)
  expect_equal(as.character(seq2),as.character(seq))
  expect_s4_class(seq2,"ModDNAString")
  #
  rnaTestSeq <- paste(alphabet(ModRNAString()), collapse = "")
  seq <- ModRNAString(rnaTestSeq)
  gr <- separate(seq)
  expect_equal(length(gr),145)
  seq2 <- combineIntoModstrings(as(seq,"RNAString"),gr)
  expect_equal(as.character(seq2),as.character(seq))
  expect_s4_class(seq2,"ModRNAString")
  #
  dnaTestSeq <- paste(alphabet(ModDNAString()), collapse = "")
  seq <- ModDNAString(dnaTestSeq)
  qseq <- PhredQuality(paste0(rep("!",
                                  length(seq)),
                              collapse = ""))
  qset <- QualityScaledModDNAStringSet(seq,
                                       qseq)
  names(qset) <- "ABC"
  gr <- separate(qset)
  actual <- combineIntoModstrings(as(qset,"DNAStringSet"),gr)
  expect_s4_class(actual,"ModDNAStringSet")
  expect_equal(as.character(actual),as.character(qset))
  actual <- combineIntoModstrings(as(qset,"DNAStringSet"),
                                gr,
                                with.qualities = TRUE)
  expect_s4_class(actual,"ModDNAStringSet")
  expect_s4_class(actual,"QualityScaledModDNAStringSet")
  expect_equal(as.character(actual),as.character(qset))
  expect_equal(unique(as.integer(quality(actual))),0)
  mcols(gr[2,])$quality <- 10
  actual <- combineIntoModstrings(as(qset,"DNAStringSet"),
                                gr,
                                with.qualities = TRUE)
  expect_s4_class(actual,"ModDNAStringSet")
  expect_s4_class(actual,"QualityScaledModDNAStringSet")
  expect_equal(as.character(actual),as.character(qset))
  expect_equal(unique(as.integer(quality(actual))),c(0,10))
  #
  rnaTestSeq <- paste(alphabet(ModRNAString()), collapse = "")
  seq <- ModRNAString(rnaTestSeq)
  qseq <- PhredQuality(paste0(rep("!",
                                  length(seq)),
                              collapse = ""))
  qset <- QualityScaledModRNAStringSet(ModRNAStringSet(list(seq,seq,seq)),
                                       c(qseq,qseq,qseq))
  names(qset) <- c("A","B","c")
  gr <- separate(qset)
  actual <- combineIntoModstrings(as(qset,"RNAStringSet"),gr)
  expect_s4_class(actual,"ModRNAStringSet")
  expect_equal(as.character(actual),as.character(qset))
  actual <- combineIntoModstrings(as(qset,"RNAStringSet"),
                                gr,
                                with.qualities = TRUE)
  expect_s4_class(actual,"ModRNAStringSet")
  expect_s4_class(actual,"QualityScaledModRNAStringSet")
  expect_equal(as.character(actual),as.character(qset))
  expect_equal(as.character(quality(actual)),as.character(quality(qset)))
  expect_equal(unique(as(quality(actual),"IntegerList")),
               IntegerList(0,0,0))
  mcols(gr[2,])$quality <- 10
  actual <- combineIntoModstrings(as(qset,"RNAStringSet"),
                                gr,
                                with.qualities = TRUE)
  expect_s4_class(actual,"ModRNAStringSet")
  expect_s4_class(actual,"QualityScaledModRNAStringSet")
  expect_equal(as.character(actual),as.character(qset))
  expect_equal(unique(as(quality(actual),"IntegerList")),
                              IntegerList(c(0,10),c(0),c(0)))
})
