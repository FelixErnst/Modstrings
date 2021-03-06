context("Sequence changes")
test_that("Sequence changes:",{
  seq <- ModDNAString("AGTC")
  expect_error(replaceLetterAt(seq,c(3,4),"CT",verbose = 1),
               "'verbose' must be TRUE or FALSE")
  expect_error(replaceLetterAt(seq,c(3,4),1),
               "'letter' must be a ModString object or a character vector")
  expect_error(replaceLetterAt(seq,c(3),"CT"),
               "total nb of letters in 'letter' must be the")
  expect_error(replaceLetterAt(seq,"c(3)","CT"),
               "'at' must be a vector of integers")
  #
  actual <- replaceLetterAt(seq,c(3,4),"CT")
  expect_equal(as.character(actual),"AGCT")
  seq <- ModRNAString("AGUC")
  actual <- replaceLetterAt(seq,c(3,4),"CU")
  expect_equal(as.character(actual),"AGCU")
  # replaceLetterAt
  expect_error(replaceLetterAt(ModDNAStringSet()),
               "'x' has no element")
  #
  set <- ModDNAStringSet(c(paste0("AGT",alphabet(ModDNAString())[12]),
                           paste0("AGT",alphabet(ModDNAString())[12])))
  expect_error(replaceLetterAt(set, 1),
               "'at' must be a matrix or list of logicals")
  expect_error(replaceLetterAt(set,
                               matrix(rep(c(FALSE,FALSE,TRUE,TRUE),1),
                                      nrow = 1,
                                      byrow = TRUE),
                               1),
               "'x' and 'at' must have the same dimensions")
  expect_error(replaceLetterAt(set,
                               matrix(rep(c(FALSE,FALSE,TRUE,TRUE),2),
                                      nrow = 2,
                                      byrow = TRUE),
                               1),
               "'letter' must be a ModStringSet object or a character vector")
  expect_error(replaceLetterAt(set,
                               matrix(rep(c(FALSE,FALSE,TRUE,TRUE),2),
                                      nrow = 2,
                                      byrow = TRUE),
                               c("CT")),
               "'x' and 'letter' must have the same length")
  expect_error(replaceLetterAt(set,
                               matrix(rep(c(0,1,2,3),2),
                                      nrow = 2,
                                      byrow = TRUE),
                               c("CT","CT")),
               "'at' must be a matrix or list of logicals")
  actual <- replaceLetterAt(set,
                            matrix(rep(c(FALSE,FALSE,TRUE,TRUE),2),
                                   nrow = 2,
                                   byrow = TRUE),
                            c("CT","CT"))
  expect_equal(as.character(actual),c("AGCT","AGCT"))
  actual <- replaceLetterAt(set,
                            list(c(FALSE,FALSE,TRUE,TRUE),
                                 c(FALSE,FALSE,TRUE,TRUE)),
                            c("CT","CT"))
  expect_equal(as.character(actual),c("AGCT","AGCT"))
  set <- ModRNAStringSet(c("AGUC","AGUC"))
  actual <- replaceLetterAt(set,
                            matrix(rep(c(FALSE,FALSE,TRUE,TRUE),2),
                                   nrow = 2,
                                   byrow = TRUE),
                            c("CU","CU"))
  expect_equal(as.character(actual),c("AGCU","AGCU"))
  # modifyNucleotides
  seq <- ModDNAString("AGTC")
  actual <- modifyNucleotides(seq,c(1,2,4),c("1mA","7mG","3mC"))
  expect_equal(as.character(actual),"\"7T'")
  seq <- ModRNAString("AGUC")
  actual <- modifyNucleotides(seq,c(1,2,3,4),c("m1A","m7G","D","m3C"))
  expect_equal(as.character(actual),"\"7D'")
  #
  expect_error(modifyNucleotides(ModDNAStringSet(),
                                 matrix(rep(c(TRUE,TRUE,FALSE,TRUE),2),
                                        nrow = 2,
                                        byrow = TRUE),
                                 list(c("1mA","7mG","3mC"),
                                      c("1mA","7mG","3mC"))),
               "'x' has no element")
  set <- ModDNAStringSet(c("AGTC","AGTC"))
  expect_error(modifyNucleotides(set,
                                 matrix(rep(c(TRUE,TRUE,FALSE,TRUE),2),
                                        nrow = 2,
                                        byrow = TRUE),
                                 list(c("1mA","7mG","3mC"))),
               "'x' and 'mod' must have the same length")
  expect_error(modifyNucleotides(set,
                                 matrix(rep(c(TRUE,TRUE,FALSE,TRUE),2),
                                        nrow = 2,
                                        byrow = TRUE),
                                 1),
               "'mod' must be a ModStringSet object, a list")
  actual <- modifyNucleotides(set,
                              matrix(rep(c(TRUE,TRUE,FALSE,TRUE),2),
                                     nrow = 2,
                                     byrow = TRUE),
                              list(c("1mA","7mG","3mC"),
                                   c("1mA","7mG","3mC")))
  expect_equal(as.character(actual),c("\"7T'","\"7T'"))
  actual <- modifyNucleotides(set,
                              list(c(TRUE,TRUE,FALSE,TRUE),
                                   c(TRUE,TRUE,FALSE,TRUE)),
                              list(c("1mA","7mG","3mC"),
                                   c("1mA","7mG","3mC")))
  expect_equal(as.character(actual),c("\"7T'","\"7T'"))
  actual <- modifyNucleotides(set,
                              matrix(rep(c(TRUE,TRUE,FALSE,TRUE),2),
                                     nrow = 2,
                                     byrow = TRUE),
                              IRanges::CharacterList(c("1mA","7mG","3mC"),
                                                     c("1mA","7mG","3mC")))
  expect_equal(as.character(actual),c("\"7T'","\"7T'"))
  actual <- modifyNucleotides(set,
                              matrix(rep(c(TRUE,TRUE,FALSE,TRUE),2),
                                     nrow = 2,
                                     byrow = TRUE),
                              ModDNAStringSet(c("\"7'","\"7'")))
  expect_equal(as.character(actual),c("\"7T'","\"7T'"))
  set <- ModRNAStringSet(c("AGUC","AGUC"))
  actual <- modifyNucleotides(set,
                              matrix(rep(c(TRUE,TRUE,TRUE,TRUE),2),
                                     nrow = 2,
                                     byrow = TRUE),
                              list(c("m1A","m7G","D","m3C"),
                                   c("m1A","m7G","D","m3C")))
  expect_equal(as.character(actual),c("\"7D'","\"7D'"))
  #
  set <- ModDNAStringSet(c("AGTC","AGTC"))
  expect_error(modifyNucleotides(set,
                              matrix(rep(c(TRUE,TRUE,FALSE,TRUE),2),
                                     nrow = 2,
                                     byrow = TRUE),
                              list(c("1mAa","7mGa","3mCa"),
                                   c("1mAa","7mGa","3mCa")),
                              nc.type = "nc"),
               "No modification for the identifiers")
  set <- ModRNAStringSet(c("AGUC","AGUC"))
  expect_error(modifyNucleotides(set,
                                 matrix(rep(c(TRUE,TRUE,TRUE,TRUE),2),
                                        nrow = 2,
                                        byrow = TRUE),
                                 list(c("1mA","7mG","dhT","3mC"),
                                      c("1mA","7mG","dhT","3mC")),
                                 nc.type = "nc"),
               "No modification for the identifiers")
  #
  expect_error(modifyNucleotides(ModDNAString("AGTC"),c(2),c("1mA")),
               "Modification type does not match the originating base")
  expect_error(modifyNucleotides(ModDNAString("AGTC"),c(2,2),c("1mA")),
               "lengths of 'at' and 'mod' need to be equal")
  expect_error(modifyNucleotides(set,
                                 matrix(rep(c(TRUE,TRUE,FALSE,TRUE),2),
                                        nrow = 2,
                                        byrow = TRUE),
                                 list(c("m1A","m7G","m3C"),
                                      c("m1A","m7G"))),
               "Dimensions of mod and 'at' must be the same")
  expect_error(modifyNucleotides(set,
                                 matrix(rep(c(TRUE,TRUE,FALSE),2),
                                        nrow = 2,
                                        byrow = TRUE),
                                 list(c("m1A","m7G","m3C"),
                                      c("m1A","m7G","m3C"))),
               "'x' and 'at' must have the same dimensions")
  #
  rnaTestSeq <- paste(alphabet(ModDNAString())[11:19], collapse = "")
  seq <- ModDNAString(rnaTestSeq)
  actual <- extractAt(seq,IRanges::IRanges(5,8))
  expect_equal(as.character(actual),"eg`b")
  #
  rnaTestSeq <- paste(alphabet(ModRNAString())[1:16], collapse = "")
  seq <- ModRNAString(rnaTestSeq)
  actual <- extractAt(seq,IRanges::IRanges(13,16))
  expect_equal(as.character(actual),"\"KO]")
  #
  dnaTestSeq <- paste(alphabet(ModDNAString())[11:19], collapse = "")
  seq <- ModDNAString(dnaTestSeq)
  actual <- replaceAt(seq,IRanges::IRanges(5,8),"AGCT")
  expect_equal(as.character(actual),"O]DJAGCTU")
  actual2 <- replaceAt(seq,c(5),"AGCT")
  expect_equal(as.character(actual2),"O]DJAGCTeg`bU")
  actual2 <- replaceAt(seq,IRanges::IRanges(5,8),DNAString("AGCT"))
  expect_equal(as.character(actual),as.character(actual2))
  actual2 <- replaceAt(seq,IRanges::IRanges(5,8),ModDNAString("AGCT"))
  expect_equal(as.character(actual),as.character(actual2))
  #
  set <- ModDNAStringSet(list(seq,seq))
  actual <- replaceAt(set,
                      List(IRanges::IRanges(5,8),IRanges::IRanges(5,8)),
                      list("AGCT","TCGA"))
  expect_equal(as.character(actual),c("O]DJAGCTU","O]DJTCGAU"))
  set <- ModDNAStringSet(list(seq,seq))
  actual2 <- replaceAt(set,
                       List(IRanges::IRanges(5,8),IRanges::IRanges(5,8)),
                       DNAStringSetList(list("AGCT","TCGA")))
  expect_equal(as.character(actual2),as.character(actual))
  #
  rnaTestSeq <- paste(alphabet(ModRNAString())[1:16], collapse = "")
  seq <- ModRNAString(rnaTestSeq)
  
  actual <- replaceAt(seq,IRanges::IRanges(7,16),"AGCUAGCUAG")
  expect_equal(as.character(actual),"ACGUN-AGCUAGCUAG")
  actual2 <- replaceAt(seq,IRanges::IRanges(7,16),"AGCTAGCTAG")
  expect_equal(as.character(actual2),"ACGUN-AGCTAGCTAG")
  actual2 <- replaceAt(seq,c(7),"AGCU")
  expect_equal(substr(as.character(actual2),1,12),"ACGUN-AGCU+.")
  actual2 <- replaceAt(seq,IRanges::IRanges(7,16),RNAString("AGCUAGCUAG"))
  expect_equal(as.character(actual2),as.character(actual))
  actual2 <- replaceAt(seq,IRanges::IRanges(7,16),ModRNAString("AGCUAGCUAG"))
  expect_equal(as.character(actual2),as.character(actual))
  set <- ModRNAStringSet(list(seq,seq))
  actual <- replaceAt(set,
                      List(IRanges::IRanges(7,16),IRanges::IRanges(7,16)),
                      list("AGCUAGCUAG","UCGAUCGAUC"))
  expect_equal(as.character(actual),c("ACGUN-AGCUAGCUAG","ACGUN-UCGAUCGAUC"))
  set <- ModRNAStringSet(list(seq,seq))
  actual2 <- replaceAt(set,
                      List(IRanges::IRanges(7,16),IRanges::IRanges(7,16)),
                      RNAStringSetList(list("AGCUAGCUAG","UCGAUCGAUC")))
  expect_equal(as.character(actual2),as.character(actual))
})
