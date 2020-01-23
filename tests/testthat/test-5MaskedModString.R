context("MaskedModString")
test_that("MaskedModString",{
  mask <- Mask(mask.width=5, start=c(2), width=c(3))
  expect_s4_class(mask,"MaskCollection")
  # MaskedModRNA
  mr <- ModRNAString("ACGU7")
  masks(mr) <- mask
  expect_s4_class(mr,"MaskedModRNAString")
  expect_equal(unmasked(mr),as(mr,"ModRNAString"))
  expect_s4_class(gaps(mr),"MaskedModRNAString")
  expect_equal(gaps(gaps(mr)),mr)
  expect_equal(as.character(mr),"A###7")
  # MaskedModDNA
  mr <- ModDNAString("ACGU7")
  masks(mr) <- mask
  expect_s4_class(mr,"MaskedModDNAString")
  expect_equal(unmasked(mr),as(mr,"ModDNAString"))
  expect_s4_class(gaps(mr),"MaskedModDNAString")
  expect_equal(gaps(gaps(mr)),mr)
  expect_equal(as.character(mr),"A###7")
})

context("ModStringViews")
test_that("ModStringViews",{
  mask <- Mask(mask.width=5, start=c(2), width=c(3))
  # ModStringsViews - RNA
  mr <- ModRNAString("ACGU7")
  masks(mr) <- mask
  expect_type(alphabetFrequency(mr),"integer")
  expect_equal(sum(alphabetFrequency(mr)),2L)
  expect_error(letterFrequency(mr),
               'argument "letters" is missing, with no default')
  expect_equal(unname(letterFrequency(mr,"7")),1L)
  v <- as(mr,"ModStringViews")
  expect_s4_class(v,"ModStringViews")
  expect_equal(consensusString(v),"?")
  expect_true(is.matrix(letterFrequency(v,"7")))
  expect_type(letterFrequency(v,"7"),"integer")
  v_D <- v
  # ModStringsViews - DNA
  mr <- ModDNAString("ACGU7")
  masks(mr) <- mask
  expect_type(alphabetFrequency(mr),"integer")
  expect_equal(sum(alphabetFrequency(mr)),2L)
  expect_error(letterFrequency(mr),
               'argument "letters" is missing, with no default')
  expect_equal(unname(letterFrequency(mr,"7")),1L)
  v <- as(mr,"ModStringViews")
  expect_s4_class(v,"ModStringViews")
  expect_equal(consensusString(v),"?")
  expect_true(is.matrix(letterFrequency(v,"7")))
  expect_type(letterFrequency(v,"7"),"integer")
  # comparison
  expect_error(v_D == v,
               "comparison between XStringViews objects")
  expect_true(all(v_D == v_D))
  expect_true(all(!v_D == ModRNAString("ACGU7")))
  expect_true(all(!v_D == RNAString("ACGU")))
})