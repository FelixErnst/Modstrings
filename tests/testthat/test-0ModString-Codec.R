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
})
