context("Matching")
test_that("Matching:",{
  dnaTestSeq <- paste(alphabet(ModDNAString()), collapse = "")
  seq <- ModDNAString(dnaTestSeq)
  actual <- matchPattern("76321",seq)
  expect_s4_class(actual,"ModStringViews")
  expect_equal(start(actual),23)
  expect_equal(end(actual),27)
  expect_equal(width(actual),5)
  actual2 <- matchPattern("7",actual)
  expect_s4_class(actual,"ModStringViews")
  expect_equal(start(actual2),23)
  expect_equal(end(actual2),23)
  expect_equal(width(actual2),1)
  #
  expect_error(matchPattern(actual,seq),
               "'pattern' must be a single string or an ModString object")
  #
  expect_equal(countPattern("76321",seq),1)
  #
  expect_equal(countPattern("76",actual),1)
  #
  set <- ModDNAStringSet(list(seq,seq))
  expect_error(matchPattern("76321",set),"please use vmatchPattern()")
  actual <- vmatchPattern("76321",set)
  expect_s4_class(actual,"MIndex")
  expect_equal(as.integer(start(actual)),c(23,23))
  expect_equal(as.integer(end(actual)),c(27,27))
  expect_equal(as.integer(width(actual)),c(5,5))
  #
  expect_equal(vcountPattern("76321",set),c(1,1))
  #
  actual <- matchLRPatterns("76321","mh",100,seq)
  expect_s4_class(actual,"ModStringViews")
  expect_equal(start(actual),23)
  expect_equal(end(actual),35)
  expect_equal(width(actual),13)
  #
  actual <- matchLRPatterns("76321","mh",5,seq)
  expect_equal(length(actual),0)
  ##############################################################################
  rnaTestSeq <- paste(alphabet(ModRNAString()), collapse = "")
  seq <- ModRNAString(rnaTestSeq)
  actual <- matchPattern("\"KO",seq)
  expect_s4_class(actual,"ModStringViews")
  expect_equal(start(actual),13)
  expect_equal(end(actual),15)
  expect_equal(width(actual),3)
  #
  actual2 <- matchPattern("\"",actual)
  expect_s4_class(actual,"ModStringViews")
  expect_equal(start(actual2),13)
  expect_equal(end(actual2),13)
  expect_equal(width(actual2),1)
  #
  expect_error(matchPattern(actual,seq),
               "'pattern' must be a single string or an ModString object")
  #
  expect_equal(countPattern("\"KO",seq),1)
  #
  expect_equal(countPattern("\"K",actual),1)
  #
  set <- ModRNAStringSet(list(seq,seq))
  expect_error(matchPattern("\"KO",set),"please use vmatchPattern()")
  actual <- vmatchPattern("\"KO",set)
  expect_s4_class(actual,"MIndex")
  expect_equal(as.integer(start(actual)),c(13,13))
  expect_equal(as.integer(end(actual)),c(15,15))
  expect_equal(as.integer(width(actual)),c(3,3))
  #
  expect_equal(vcountPattern("\"KO",set),c(1,1))
  #
  actual <- matchLRPatterns("\"KO","RL",100,seq)
  expect_s4_class(actual,"ModStringViews")
  expect_equal(start(actual),13)
  expect_equal(end(actual),105)
  expect_equal(width(actual),93)
  #
  actual <- matchLRPatterns("\"KO","5F",50,seq)
  expect_s4_class(actual,"ModStringViews")
  expect_equal(length(actual),0)
  #
  actual <- matchLRPatterns("\"KO","5F",200,seq)
  expect_s4_class(actual,"ModStringViews")
  expect_equal(start(actual),13)
  expect_equal(end(actual),82)
  expect_equal(width(actual),70)
})
