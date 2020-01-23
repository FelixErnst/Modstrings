context("MaskedModString")
test_that("MaskedModString",{
  mask <- Mask(mask.width=5, start=c(2), width=c(3))
  expect_s4_class(mask,"MaskCollection")
  # MaskedModRNA
  mr <- ModRNAString("ACGU7")
  masks(mr) <- mask
  expect_s4_class(mr,"MaskedModRNAString")
  expect_equal(unmasked(mr),as(mr,"ModRNAString"))
  v <- as(mr,"ModStringViews")
  expect_s4_class(v,"ModStringViews")
  expect_s4_class(gaps(mr),"MaskedModRNAString")
  expect_equal(gaps(gaps(mr)),mr)
  expect_equal(as.character(mr),"A###7")
  # MaskedModDNA
  mr <- ModDNAString("ACGU7")
  masks(mr) <- mask
  expect_s4_class(mr,"MaskedModDNAString")
  expect_equal(unmasked(mr),as(mr,"ModDNAString"))
  v <- as(mr,"ModStringViews")
  expect_s4_class(v,"ModStringViews")
  expect_s4_class(gaps(mr),"MaskedModDNAString")
  expect_equal(gaps(gaps(mr)),mr)
  expect_equal(as.character(mr),"A###7")
})