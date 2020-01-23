context("Writing and Reading ModStringSet")
test_that("Writing and Reading ModStringSet:",{
  seqs <- paste0(paste(alphabet(ModDNAString()), collapse = ""),
                 c("A","G","C"))
  set <- ModDNAStringSet(seqs)
  names(set) <- c("A","B","C")
  file <- tempfile()
  write("",file = file)
  expect_error(readModDNAStringSet(file),
               "No fasta blocks found.")
  expect_error(writeModStringSet(),
               "'x' must be an ModStringSet object")
  expect_error(writeModStringSet(set),
               "'filepath' must be a single string")
  expect_error(writeModStringSet(set, file, format = 1),
               "'format' must be a single string")
  expect_error(writeModStringSet(set, file, append = 1),
               "'append' must be a single logical")
  writeModStringSet(set, file)
  read <- readModDNAStringSet(file)
  expect_equal(names(readModDNAStringSet(file,nrec = 2)),c("A","B"))
  expect_equal(names(readModDNAStringSet(file,nrec = 2, skip = 1)),c("B","C"))
  expect_equal(names(readModDNAStringSet(file,nrec = 2, skip = 2)),c("C"))
  # expect_equal does not work because of encoding
  expect_true(all(as.character(set) == as.character(read)))
  expect_equal(names(set),names(read))
  expect_error(readModDNAStringSet(file, format = "fastq"),
               "Wrong format. '@' not found. Aborting")
  expect_error(readModDNAStringSet(file, with.qualities = TRUE,
                                   format = "fasta"),
               "The 'with.qualities' argument is only supported when")
  writeModStringSet(set, file, format = "fastq")
  expect_error(readModDNAStringSet(file, format = "fasta"),
               "reading FASTA file")
  read <- readModDNAStringSet(file, format = "fastq")
  # expect_equal does not work because of encoding
  expect_true(all(as.character(set) == as.character(read)))
  expect_equal(names(set),names(read))
  expect_equal(readModDNAStringSet(file, format = "fastq",
                                   seek.first.rec = TRUE),
               read)
  expect_error(readModDNAStringSet(file, format = "fastq",
                                   seek.first.rec = 1),
               "'seek.first.rec' must be TRUE or FALSE")
  expect_error(readModDNAStringSet(file, format = "fastq",
                                   with.qualities = 1),
               "'with.qualities' must be TRUE or FALSE")
  actual <- readModDNAStringSet(file, format = "fastq",
                                with.qualities = TRUE)
  expect_s4_class(actual,"ModDNAStringSet")
  expect_s4_class(mcols(actual),"DataFrame")
  expect_s4_class(mcols(actual)$qualities,"BStringSet")
  qualities <- mcols(actual)$qualities
  actual <- readModDNAStringSet(file, format = "fastq",
                                with.qualities = TRUE, nrec = 2)
  expect_equal(names(actual),c("A","B"))
  actual <- readModDNAStringSet(file, format = "fastq",
                                with.qualities = TRUE, nrec = 2, skip = 1)
  expect_equal(names(actual),c("B","C"))
  actual <- readModDNAStringSet(file, format = "fastq",
                                with.qualities = TRUE, nrec = 2, skip = 2)
  expect_equal(names(actual),c("C"))
  expect_error(writeModStringSet(set, file, format = "fastq", qualities = 1),
                   "'qualities' must be NULL or a BStringSet object")
  expect_error(writeModStringSet(set, file, format = "fastq",
                                     qualities = qualities[1]),
                   "'x' and 'qualities' must have the same length")
  expect_invisible(writeModStringSet(set, file, format = "fastq"))
  read <- readQualityScaledModDNAStringSet(file)
  expect_s4_class(read,"QualityScaledModDNAStringSet")
  writeModStringSet(set, file, format = "fastq", qualities = qualities)
  read2 <- readQualityScaledModDNAStringSet(file)
  expect_equal(read,read2)
  expect_error(writeQualityScaledModStringSet(set, file),
               "'x' must be a QualityScaledXStringSet object")
  writeQualityScaledModStringSet(read, file)
  read2 <- readQualityScaledModDNAStringSet(file)
  expect_equal(read,read2)
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
