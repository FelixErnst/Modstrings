context("Codec")
test_that("Codec:",{
  dnacodec <- Modstrings:::MOD_DNA_BASE_CODES
  rnacodec <- Modstrings:::MOD_RNA_BASE_CODES
  dnacodec <- dnacodec[dnacodec$abbrev != "",]
  rnacodec <- rnacodec[rnacodec$abbrev != "",]
  expect_true(!any(duplicated(dnacodec$name)))
  expect_true(!any(duplicated(dnacodec$short_name)))
  expect_true(!any(duplicated(dnacodec$nc)))
  expect_true(!any(duplicated(dnacodec$abbrev)))
  expect_true(!any(duplicated(dnacodec$value)))
  expect_true(!any(duplicated(dnacodec$oneByteLetter)))
  expect_true(!any(duplicated(rnacodec$name)))
  expect_true(!any(duplicated(rnacodec$short_name)))
  expect_true(!any(duplicated(rnacodec$nc)))
  expect_true(!any(duplicated(rnacodec$abbrev)))
  expect_true(!any(duplicated(rnacodec$value)))
  expect_true(!any(duplicated(rnacodec$oneByteLetter)))
  actual <- Modstrings:::.load_mod_dictionary("test_DNA_codes.txt")
  expect_s4_class(actual,"DataFrame")
  expect_equal(colnames(actual),c("name","short_name","abbrev","orig_base",
                                  "nc","value","oneByteLetter"))
  #
  actual <- Modstrings:::.new_ModStringCodec(Modstrings:::MOD_DNA_BASE_CODES,
                                             c(Modstrings:::.DNA_BASE_CODES, 
                                               Modstrings:::additional_base_codes))
  expect_s4_class(actual,"ModStringCodec")
  expect_type(Modstrings:::letters(actual),"character")
  expect_type(Modstrings:::oneByteCodes(actual),"character")
  expect_type(Modstrings:::conversion(actual),"logical")
  expect_type(Modstrings:::originatingBase(actual),"character")
  expect_equal(unique(Modstrings:::originatingBase(actual)),c("A","C","T","G",
                                                              "N","-","+","."))
  expect_type(Modstrings:::values(actual),"integer")
  expect_type(Modstrings:::lettersEscaped(actual),"character")
  expect_type(Modstrings:::oneByteCodesEscaped(actual),"character")
  expect_type(Modstrings:::lettersNeedEscape(actual),"logical")
  expect_type(Modstrings:::oneByteCodesNeedEscape(actual),"logical")
  expect_s4_class(Modstrings:::additionalInfo(actual),"DataFrame")
})
