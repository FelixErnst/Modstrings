context("Writing and Reading ModStringSet")
test_that("Writing and Reading ModStringSet:",{
  seqs <- paste0(paste(alphabet(ModDNAString()), collapse = ""),
                 c("A","G","C"))
  set <- ModDNAStringSet(seqs)
  names(set) <- c("A","B","C")
  file <- tempfile()
  writeModStringSet(set, file)
  read <- readModDNAStringSet(file)
  # expect_equal does not work because of encoding
  expect_true(all(as.character(set) == as.character(read)))
  expect_equal(names(set),names(read))
  writeModStringSet(set, file, format = "fastq")
  read <- readModDNAStringSet(file, format = "fastq")
  # expect_equal does not work because of encoding
  expect_true(all(as.character(set) == as.character(read)))
  expect_equal(names(set),names(read))
  read <- readQualityScaledModDNAStringSet(file)
  expect_s4_class(read,"QualityScaledModDNAStringSet")
  ####################
  seqs <- paste0(paste(alphabet(ModRNAString()), collapse = ""),
                 c("A","G","C"))
  set <- ModRNAStringSet(seqs)
  names(set) <- c("A","B","C")
  file <- tempfile()
  writeModStringSet(set, file)
  read <- readModRNAStringSet(file)
  # expect_equal does not work because of encoding
  expect_true(all(as.character(set) == as.character(read)))
  expect_equal(names(set),names(read))
  writeModStringSet(set, file, format = "fastq")
  read <- readModRNAStringSet(file, format = "fastq")
  # expect_equal does not work because of encoding
  expect_true(all(as.character(set) == as.character(read)))
  expect_equal(names(set),names(read))
  read <- readQualityScaledModRNAStringSet(file)
  expect_s4_class(read,"QualityScaledModRNAStringSet")
})
