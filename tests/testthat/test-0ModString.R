context("ModString")
test_that("ModString:",{
  expect_equal(as.character(ModDNAString("AGCT")),"AGCT")
  expect_equal(as.character(ModRNAString("AGCT")),"AGCT")
  
  dnaTestSeq <- paste(alphabet(ModDNAString()), collapse = "")
  expect_equal(as.character(ModDNAString(dnaTestSeq)),dnaTestSeq)
  expect_equal(length(ModDNAString(dnaTestSeq)),nchar(dnaTestSeq))
  expect_output(show(ModDNAString(dnaTestSeq)))
  dnaTestSeq <- paste0(dnaTestSeq,dnaTestSeq,dnaTestSeq)
  expect_equal(as.character(ModDNAString(dnaTestSeq)),dnaTestSeq)
  rnaTestSeq <- paste(alphabet(ModRNAString()), collapse = "")
  expect_equal(as.character(ModRNAString(rnaTestSeq)),rnaTestSeq)
  expect_output(show(ModRNAString(rnaTestSeq)))
  expect_equal(length(ModRNAString(rnaTestSeq)),nchar(rnaTestSeq))
  rnaTestSeq <- paste0(rnaTestSeq,rnaTestSeq,rnaTestSeq)
  expect_equal(as.character(ModRNAString(rnaTestSeq)),rnaTestSeq)
})
context("ModStringSet")
test_that("ModStringSet:",{
  seqs <- paste0(paste(alphabet(ModDNAString()), collapse = ""),
                 c("A","G","C"))
  names(seqs) <- c("A","B","C")
  actual <- ModDNAStringSet(seqs)
  expect_output(show(actual))
  expect_named(actual,names(seqs))
  expect_equal(as.character(actual[1]),seqs[1])
  expect_equal(as.character(actual[2]),seqs[2])
  expect_equal(as.character(actual[3]),seqs[3])
  expect_equal(as.character(ModDNAStringSet(list(ModDNAString(seqs[1]),
                                                 ModDNAString(seqs[2]))[2])),
               unname(seqs[2]))
  
  seqs <- paste0(paste(alphabet(ModRNAString()), collapse = ""),
                 c("A","G","C"))
  names(seqs) <- c("A","B","C")
  actual <- ModRNAStringSet(seqs)
  expect_output(show(actual))
  expect_named(actual,names(seqs))
  expect_equal(as.character(actual[1]),seqs[1])
  expect_equal(as.character(actual[2]),seqs[2])
  expect_equal(as.character(actual[3]),seqs[3])
  expect_equal(as.character(ModRNAStringSet(list(ModRNAString(seqs[1]),
                                                 ModRNAString(seqs[2]))[2])),
               unname(seqs[2]))
})
context("ModStringSetList")
test_that("ModStringSetList:",{
  seqs1 <- paste0(paste(alphabet(ModDNAString()), collapse = ""),
                  c("A","G","C"))
  seqs2 <- paste0(paste(alphabet(ModDNAString()), collapse = ""),
                  c("C","T","A"))
  actual <- ModDNAStringSetList(list(ModDNAStringSet(seqs1),
                                     ModDNAStringSet(seqs2)))
  expect_output(show(actual))
  expect_equal(as.character(actual[[1]][1]),seqs1[1])
  expect_equal(as.character(actual[[2]][1]),seqs2[1])
  expect_equal(as.character(actual[[1]][3]),seqs1[3])
  expect_equal(as.character(actual[[2]][2]),seqs2[2])
  seqs1 <- paste0(paste(alphabet(ModRNAString()), collapse = ""),
                  c("A","G","C"))
  seqs2 <- paste0(paste(alphabet(ModRNAString()), collapse = ""),
                  c("C","T","A"))
  actual <- ModRNAStringSetList(list(ModRNAStringSet(seqs1),
                                     ModRNAStringSet(seqs2)))
  expect_output(show(actual))
  expect_equal(as.character(actual[[1]][1]),seqs1[1])
  expect_equal(as.character(actual[[2]][1]),seqs2[1])
  expect_equal(as.character(actual[[1]][3]),seqs1[3])
  expect_equal(as.character(actual[[2]][2]),seqs2[2])
})

context("ModStringViews")
test_that("ModStringViews:",{
  seq <- ModDNAString(paste(alphabet(ModDNAString()),
                            collapse = ""))
  v <- Views(seq, start = 3:1, end = 26:28)
  expect_output(show(v))
  expect_s4_class(v, "ModStringViews")
  expect_type(as.character(v), "character")
  expect_equal(length(v), 3)
  expect_equal(start(v), c(3,2,1))
  expect_equal(end(v), c(26,27,28))
  expect_equal(width(v), c(24,26,28))
  set <- as(v,"XStringSet")
  expect_s4_class(set, "ModDNAStringSet")
  expect_equal(length(set), 3)
  expect_equal(width(set), c(24,26,28))
  #
  seq <- ModRNAString(paste(alphabet(ModRNAString()),
                            collapse = ""))
  v <- Views(seq, start = 3:1, end = 26:28)
  expect_output(show(v))
  expect_s4_class(v, "ModStringViews")
  expect_type(as.character(v), "character")
  expect_equal(length(v), 3)
  expect_equal(start(v), c(3,2,1))
  expect_equal(end(v), c(26,27,28))
  expect_equal(width(v), c(24,26,28))
  set <- as(v,"XStringSet")
  expect_s4_class(set, "ModRNAStringSet")
  expect_equal(length(set), 3)
  expect_equal(width(set), c(24,26,28))
})


context("Conversion to Biostrings")
test_that("Conversion to Biostrings:",{
  ##############################################################################
  # ModString
  seq <- ModDNAString("AGCT")
  dna <- DNAString(seq)
  expect_equal(as.character(seq),as.character(dna))
  expect_s4_class(as(dna,"ModDNAString"),"ModDNAString")
  seqtype(seq) <- "DNA"
  expect_equal(seq,dna)
  #
  seq <- ModRNAString("AGCU")
  rna <- RNAString(seq)
  expect_equal(as.character(seq),as.character(rna))
  expect_s4_class(as(rna,"ModRNAString"),"ModRNAString")
  seqtype(seq) <- "RNA"
  expect_equal(seq,rna)
  #
  seq <- ModRNAString("AGCT")
  rna <- RNAString(seq)
  expect_equal(as.character(rna),"AGCU")
  #
  seq <- paste(alphabet(ModDNAString()), collapse = "")
  seq <- ModDNAString(seq)
  dna <- DNAString(seq)
  expect_equal(as.character(dna),
               paste0("ACGTN-+.TTTTTTTTTTTTTAGGGGGGGGGGGCCCCCCCCCCCAAAAAAAAAA",
                      "A"))
  expect_equal(as.character(dna),
               as.character(ModDNAString(dna)))
  #
  seq <- paste(alphabet(ModRNAString()), collapse = "")
  seq <- ModRNAString(seq)
  rna <- RNAString(seq)
  expect_equal(
    as.character(rna),
    paste0("ACGUN-+.AGAUAGAUAUCAAAAAAAAUUCUACGAUUUAGUUUUCUUGUCUUUUUUUUUUUUUUUU",
           "UUUUUUUCCCCUUUUUUUUUUCUUUUGGGGGGAGGGGGGGCCCCCCAAAAAAAAAAAAACGAUGG",
           "GAGAGGGGGUGGNNACGUUUGG"))
  expect_equal(as.character(rna),
               as.character(ModRNAString(rna)))
  ##############################################################################
  # ModStringSet
  seqs <- paste0(paste(alphabet(ModDNAString()), collapse = ""),
                 c("A","G","C"))
  set <- ModDNAStringSet(seqs)
  names(set) <- c("A","B","C")
  dnaSeqs <- c(
    A = paste0("ACGTN-+.TTTTTTTTTTTTTAGGGGGGGGGGGCCCCCCCCCCCAAAAAAAAAA",
               "AA"),
    B = paste0("ACGTN-+.TTTTTTTTTTTTTAGGGGGGGGGGGCCCCCCCCCCCAAAAAAAAAA",
               "AG"),
    C = paste0("ACGTN-+.TTTTTTTTTTTTTAGGGGGGGGGGGCCCCCCCCCCCAAAAAAAAAA",
               "AC"))
  expect_equal(as.character(DNAStringSet(set)),
               dnaSeqs)
  expect_equal(as.character(DNAStringSet(set)),
               as.character(ModDNAStringSet(DNAStringSet(set))))
  #
  seqs <- paste0(paste(alphabet(ModRNAString()), collapse = ""),
                 c("A","G","C"))
  set <- ModRNAStringSet(seqs)
  names(set) <- c("A","B","C")
  rnaSeqs <- c(
    A = paste0("ACGUN-+.AGAUAGAUAUCAAAAAAAAUUCUACGAUUUAGUUUUCUUGUCUUUUUUUUUUUU",
               "UUUUUUUUUUUCCCCUUUUUUUUUUCUUUUGGGGGGAGGGGGGGCCCCCCAAAAAAAAAAAA",
               "ACGAUGGGAGAGGGGGUGGNNACGUUUGGA"),
    B = paste0("ACGUN-+.AGAUAGAUAUCAAAAAAAAUUCUACGAUUUAGUUUUCUUGUCUUUUUUUUUUUU",
               "UUUUUUUUUUUCCCCUUUUUUUUUUCUUUUGGGGGGAGGGGGGGCCCCCCAAAAAAAAAAAA",
               "ACGAUGGGAGAGGGGGUGGNNACGUUUGGG"),
    C = paste0("ACGUN-+.AGAUAGAUAUCAAAAAAAAUUCUACGAUUUAGUUUUCUUGUCUUUUUUUUUUUU",
               "UUUUUUUUUUUCCCCUUUUUUUUUUCUUUUGGGGGGAGGGGGGGCCCCCCAAAAAAAAAAAA",
               "ACGAUGGGAGAGGGGGUGGNNACGUUUGGC"))
  expect_equal(as.character(RNAStringSet(set)),
               rnaSeqs)
  expect_equal(as.character(RNAStringSet(set)),
               as.character(ModRNAStringSet(RNAStringSet(set))))
  #
  seq <- ModDNAString(paste(alphabet(ModDNAString()),
                            collapse = ""))
  v <- Views(seq, start = 3:1, end = 6:8)
  expect_equal(as.character(DNAStringSet(v))[1],"GTN-")
  seq <- ModRNAString(paste(alphabet(ModRNAString()),
                            collapse = ""))
  v <- Views(seq, start = 3:1, end = 6:8)
  expect_equal(as.character(RNAStringSet(v))[1],"GUN-")
  #
  seq <- ModDNAString("AGCT")
  expect_equal(as.character(ModRNAString(seq)),"AGCU")
  seq <- ModRNAString("AGCT")
  expect_equal(as.character(ModDNAString(seq)),"AGCT")
  seq <- ModRNAString("AGCU")
  expect_equal(as.character(ModDNAString(seq)),"AGCT")
  #
  seqs <- paste0(c("AGCT","AGCT","AGCT"),
                 c("A","G","C"))
  set <- ModDNAStringSet(seqs)
  expect_equal(as.character(ModRNAStringSet(set))[1],"AGCUA")
  seqs <- paste0(c("AGCU","AGCU","AGCU"),
                 c("A","G","C"))
  set <- ModRNAStringSet(seqs)
  expect_equal(as.character(ModDNAStringSet(set))[1],"AGCTA")
  seqs <- paste0(c("AGCT","AGCT","AGCT"),
                 c("A","G","C"))
  set <- ModRNAStringSet(seqs)
  expect_equal(as.character(ModDNAStringSet(set))[1],"AGCTA")
  #
  seq <- ModDNAString("AGCT")
  expect_equal(as.character(ModDNAStringSet(seq)),"AGCT")
  expect_equal(as.character(ModRNAStringSet(seq)),"AGCU")
  seq <- ModRNAString("AGCU")
  expect_equal(as.character(ModDNAStringSet(seq)),"AGCT")
  expect_equal(as.character(ModRNAStringSet(seq)),"AGCU")
  ##############################################################################
  # errors
  expect_output(expect_error(ModDNAString(paste(alphabet(ModRNAString()),
                                                collapse = "")),
                             "Invalid char"))
  expect_output(expect_error(ModDNAString(paste(alphabet(RNAString()),
                                                collapse = "")),
                             "Invalid char"))
})
