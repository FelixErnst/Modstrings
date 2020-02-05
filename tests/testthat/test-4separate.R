context("ModString separate/combine")
test_that("ModString separate/combine:",{
  dnaTestSeq <- paste(alphabet(ModDNAString()), collapse = "")
  seq <- ModDNAString(dnaTestSeq)
  gr <- separate(seq)
  #
  gr2 <- gr
  strand(gr2[1]) <- "-"
  expect_warning(combineIntoModstrings(as(seq,"DNAString"),gr2),
                 "Annotations on the minus strand will not be used.")
  gr2 <- gr
  end(gr2[2]) <- 60L 
  start(gr2[2]) <- 60L 
  expect_error(expect_warning(combineIntoModstrings(as(seq,"DNAString"),gr2)),
               "GRanges object contains coordinates out of bounds")
  start(gr2[2]) <- 10L
  end(gr2[2]) <- 11L
  expect_error(combineIntoModstrings(as(seq,"DNAString"),gr2),
               "width\\(\\) of GRanges elements must all be == 1")
  gr2 <- gr
  mcols(gr2)$mod <- NULL
  expect_error(combineIntoModstrings(as(seq,"DNAString"),gr2),
               "GRanges object does not contain a 'mod' column.")
  #
  expect_equal(length(gr),47)
  seq2 <- combineIntoModstrings(as(seq,"DNAString"),gr)
  expect_equal(as.character(seq2),as.character(seq))
  expect_s4_class(seq2,"ModDNAString")
  #
  rnaTestSeq <- paste(alphabet(ModRNAString()), collapse = "")
  seq <- ModRNAString(rnaTestSeq)
  gr <- separate(seq)
  #
  gr2 <- gr
  gr2 <- gr2[c(5,24)]
  start(gr2) <- 13L
  end(gr2) <- 13L
  #
  expect_error(combineIntoModstrings(as(seq,"RNAString"),shift(gr2,1)),
               "Modification type does not match")
  expect_warning(expect_equal(as.character(combineIntoModstrings(as(seq,"RNAString"),shift(gr2,1),
                                                                 stop.on.error = FALSE)),
                              as.character(as(seq,"RNAString"))),
                 "Modification type does not match")
  #
  actual <- combineIntoModstrings(as(seq,"RNAString"),gr2)
  expect_equal(length(separate(actual)),1L)
  expect_equal(unname(separate(actual)$mod),"m1Am")
  mcols(gr2)$quality <- 10L
  expect_error(combineIntoModstrings(as(seq,"RNAString"),gr2),
               "Multiple modifications found for position '13'.")
  #
  expect_equal(length(gr),145)
  seq2 <- combineIntoModstrings(as(seq,"RNAString"),gr)
  expect_equal(as.character(seq2),as.character(seq))
  expect_s4_class(seq2,"ModRNAString")
  #
  dnaTestSeq <- paste(alphabet(ModDNAString()), collapse = "")
  set <- ModDNAStringSet(c("A" = dnaTestSeq,
                           "B" = dnaTestSeq,
                           "C" = dnaTestSeq))
  gr <- separate(set)
  expect_equal(length(gr),141)
  set2 <- combineIntoModstrings(as(set,"DNAStringSet"),gr)
  expect_equal(as.character(set2),as.character(set))
  expect_s4_class(set2,"ModDNAStringSet")
  grl <- split(gr,seqnames(gr))
  set3 <- combineIntoModstrings(as(set,"DNAStringSet"),grl)
  expect_equal(as.character(set3),as.character(set))
  expect_s4_class(set3,"ModDNAStringSet")
  set2 <- combineIntoModstrings(as(set,"DNAStringSet"),gr[seqnames(gr) %in% c("A","B")])
  expect_s4_class(set2,"ModDNAStringSet")
  expect_true(as(set,"DNAStringSet")[3] == set2[3])
  #
  rnaTestSeq <- paste(alphabet(ModRNAString()), collapse = "")
  set <- ModRNAStringSet(c("A" = rnaTestSeq,
                           "B" = rnaTestSeq,
                           "C" = rnaTestSeq))
  gr <- separate(set)
  expect_equal(length(gr),435)
  #
  expect_error(combineIntoModstrings(as(set,"RNAStringSet"),
                                     GenomicRanges::GRangesList()),
               "'gr' is empty \\(length\\(\\) == 0L\\)")
  expect_error(combineIntoModstrings(unname(as(set,"RNAStringSet")),gr),
               "Names of XStringSet object do not match the seqnames")
  gr2 <- gr
  end(gr2[1]) <- 200L
  start(gr2[1]) <- 200L
  expect_error( combineIntoModstrings(as(set,"RNAStringSet"),gr2),
                "Elements of the GRangesList object contain coordinates out")
  gr2 <- gr
  end(gr2[1]) <- 10L
  expect_error(combineIntoModstrings(as(set,"RNAStringSet"),gr2),
               "width\\(\\) of GRangesList elements must all be == 1")
  gr2 <- gr
  gr2 <- gr2[c(5,24)]
  start(gr2) <- 13L
  end(gr2) <- 13L
  actual <- combineIntoModstrings(as(set,"RNAStringSet"),gr2)
  expect_equal(length(separate(actual)),1L)
  expect_equal(unname(separate(actual)$mod),"m1Am")
  mcols(gr2)$quality <- 10L
  expect_error(combineIntoModstrings(as(seq,"RNAString"),gr2),
               "Multiple modifications found for position '13'.")
  #
  set2 <- combineIntoModstrings(as(set,"RNAStringSet"),gr)
  expect_equal(as.character(set2),as.character(set))
  expect_s4_class(set2,"ModRNAStringSet")
  set2 <- combineIntoModstrings(as(set,"RNAStringSet"),
                                gr[seqnames(gr) %in% c("A","B")])
  expect_s4_class(set2,"ModRNAStringSet")
  expect_true(as(set,"RNAStringSet")[3] == set2[3])
  grl <- split(gr,seqnames(gr))
  set3 <- combineIntoModstrings(as(set,"RNAStringSet"),grl)
  expect_equal(as.character(set3),as.character(set))
  expect_s4_class(set3,"ModRNAStringSet")
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
  names(qset) <- c("A","B","C")
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
  actual <- combineIntoModstrings(as(qset,"RNAStringSet"),
                                  gr[seqnames(gr) %in% c("A","B")],
                                  with.qualities = TRUE)
  expect_s4_class(actual,"ModRNAStringSet")
  expect_s4_class(actual,"QualityScaledModRNAStringSet")
})
