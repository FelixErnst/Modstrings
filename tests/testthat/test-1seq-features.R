context("Sequence features")
test_that("hasOnlyBaseLetters:",{
  seq <- paste(alphabet(ModDNAString()), collapse = "")
  seq2 <- ModDNAString(seq)
  expect_equal(letter(seq,10), letter(seq2,10))
  seq <- paste(alphabet(ModRNAString()), collapse = "")
  seq2 <- ModRNAString(seq)
  expect_equal(letter(seq,10), letter(seq2,10))
})
test_that("hasOnlyBaseLetters:",{
  seq <- paste(alphabet(ModDNAString()), collapse = "")
  seq <- ModDNAString(seq)
  expect_false(hasOnlyBaseLetters(seq))
  seq <- "AGCT"
  seq <- ModDNAString(seq)
  expect_true(hasOnlyBaseLetters(seq))
  #
  seq <- paste(alphabet(ModRNAString()), collapse = "")
  seq <- ModRNAString(seq)
  expect_false(hasOnlyBaseLetters(seq))
  seq <- "AGCU"
  seq <- ModRNAString(seq)
  expect_true(hasOnlyBaseLetters(seq))
  seq <- "AGCT"
  seq <- ModRNAString(seq)
  expect_true(hasOnlyBaseLetters(seq))
  #
  seq <- paste(alphabet(ModDNAString()), collapse = "")
  set <- ModDNAStringSet(seq)
  expect_false(hasOnlyBaseLetters(set))
  seq <- "AGCT"
  set <- ModDNAStringSet(seq)
  expect_true(hasOnlyBaseLetters(set))
  #
  seq <- paste(alphabet(ModRNAString()), collapse = "")
  set <- ModRNAStringSet(seq)
  expect_false(hasOnlyBaseLetters(set))
  seq <- "AGCU"
  set <- ModRNAStringSet(seq)
  expect_true(hasOnlyBaseLetters(set))
  seq <- "AGCT"
  set <- ModRNAStringSet(seq)
  expect_true(hasOnlyBaseLetters(set))
})

test_that("alphabetFrequency:",{
  seq <- paste(alphabet(ModDNAString()), collapse = "")
  seq <- ModDNAString(seq)
  actual <- alphabetFrequency(seq)
  expect_equal(sum(actual),length(seq))
  #
  seq <- paste(alphabet(ModRNAString()), collapse = "")
  seq <- ModRNAString(seq)
  actual <- alphabetFrequency(seq)
  expect_equal(sum(actual),length(seq))
  #
  seq <- paste(alphabet(ModDNAString()), collapse = "")
  set <- ModDNAStringSet(c(seq,seq))
  actual <- alphabetFrequency(set)
  expect_equal(nrow(actual),length(set))
  expect_equal(ncol(actual),width(set)[1])
  #
  seq <- paste(alphabet(ModRNAString()), collapse = "")
  set <- ModRNAStringSet(c(seq,seq))
  actual <- alphabetFrequency(set)
  expect_equal(nrow(actual),length(set))
  expect_equal(ncol(actual),width(set)[1])
})

test_that("letterFrequency:",{
  seq <- paste(alphabet(ModDNAString()), collapse = "")
  seq <- ModDNAString(seq)
  actual <- letterFrequencyInSlidingView(seq,
                                         10,
                                         "AC")
  expect_equal(sum(actual),3)
  expect_equal(nrow(actual),length(seq) - 9)
  #
  seq <- paste(alphabet(ModRNAString()), collapse = "")
  seq <- ModRNAString(seq)
  actual <- letterFrequencyInSlidingView(seq,
                                         10,
                                         "AC")
  expect_equal(sum(actual),3)
  expect_equal(nrow(actual),length(seq) - 9)
  #
  seq <- paste(alphabet(ModDNAString()), collapse = "")
  set <- ModDNAStringSet(c(seq,seq))
  actual <- letterFrequency(set,
                            "AC")
  expect_equal(nrow(actual),length(set))
  expect_equal(sum(actual),4)
  #
  seq <- paste(alphabet(ModRNAString()), collapse = "")
  set <- ModRNAStringSet(c(seq,seq))
  actual <- letterFrequency(set,
                            "AC")
  expect_equal(nrow(actual),length(set))
  expect_equal(sum(actual),4)
})

test_that("consensusMatrix/consensusString:",{
  seq <- paste(alphabet(ModDNAString()), collapse = "")
  set <- ModDNAStringSet(c(seq,seq))
  actual <- consensusMatrix(set)
  expect_equal(nrow(actual),width(set)[1])
  expect_equal(ncol(actual),width(set)[1])
  actual <- consensusString(set)
  expect_equal(actual,seq)
  #
  seq <- paste(alphabet(ModRNAString()), collapse = "")
  set <- ModRNAStringSet(c(seq,seq))
  actual <- consensusMatrix(set)
  expect_equal(nrow(actual),width(set)[1])
  expect_equal(ncol(actual),width(set)[1])
  actual <- consensusString(set)
  expect_equal(actual,seq)
})

test_that("length/width/narrow/etc:",{
  seq <- ModDNAString(paste(alphabet(ModDNAString()), collapse = ""))
  l <- length(alphabet(ModDNAString()))
  expect_equal(length(seq),l)
  set <- ModDNAStringSet(paste(alphabet(ModDNAString()), collapse = ""))
  expect_equal(length(set),1)
  expect_equal(width(set),l)
  set <- c(set,set)
  expect_equal(width(set),c(l,l))
  #
  seq <- ModRNAString(paste(alphabet(ModRNAString()), collapse = ""))
  l <- length(alphabet(ModRNAString()))
  expect_equal(length(seq),l)
  set <- ModRNAStringSet(paste(alphabet(ModRNAString()), collapse = ""))
  expect_equal(length(set),1)
  expect_equal(width(set),l)
  set <- c(set,set)
  expect_equal(width(set),c(l,l))
  #
  seq <- ModDNAString(paste(alphabet(ModDNAString()), collapse = ""))
  seq2 <- subseq(seq,1,20)
  seq3 <- ModDNAString(paste(alphabet(ModDNAString())[seq_len(20)],
                             collapse = ""))
  expect_equal(as.character(seq2),as.character(seq3))
  #
  set <- ModDNAStringSet(paste(alphabet(ModDNAString()), collapse = ""))
  set2 <- subseq(set,1,20)
  set3 <- ModDNAStringSet(paste(alphabet(ModDNAString())[seq_len(20)],
                                collapse = ""))
  expect_equal(as.character(set2),as.character(set3))
  #
  seq <- ModRNAString(paste(alphabet(ModRNAString()), collapse = ""))
  seq2 <- subseq(seq,1,20)
  seq3 <- ModRNAString(paste(alphabet(ModRNAString())[seq_len(20)],
                             collapse = ""))
  expect_equal(as.character(seq2),as.character(seq3))
  #
  set <- ModRNAStringSet(paste(alphabet(ModRNAString()), collapse = ""))
  set2 <- subseq(set,1,20)
  set3 <- ModRNAStringSet(paste(alphabet(ModRNAString())[seq_len(20)],
                                collapse = ""))
  expect_equal(as.character(set2),as.character(set3))
  
})
