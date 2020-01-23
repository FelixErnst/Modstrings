context("MaskedModString")
test_that("MaskedModString",{
  chr <- "AGC@"
  expect_output(expect_error(ModRNAString(chr),
                             "Invalid character\\(s\\) - see above"))
  expect_type(sanitizeFromModomics(chr),"character")
  expect_equal(sanitizeFromModomics(chr),"AGC÷")
  #
  chr <- "AGC+"
  expect_s4_class(ModRNAString(chr),"ModRNAString")
  expect_type(sanitizeFromtRNAdb(chr),"character")
})
