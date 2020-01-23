context("Sequence comparison")
test_that("Sequence comparison:",{
  r <- DNAString("ACGTG")
  mr <- modifyNucleotides(r,5,"7mG")
  #
  r <- RNAString("ACGUG")
  rs <- RNAStringSet(list(r,r,r,r,r))
  names(rs) <- paste0("Sequence", seq_along(rs))
  mr <- modifyNucleotides(r,5,"m7G")
  expect_true(r == ModRNAString(r))
  expect_true(r == RNAString(mr))
  mrs <- ModRNAStringSet(rs)
  expect_true(all(rs == mrs))
  gr <- GRanges(seqnames = names(rs)[c(1,1,2,3,3,4,5,5)],
                 ranges = IRanges(start = c(4,5,5,4,5,5,4,5),width = 1),
                 mod = c("D","m7G","m7G","D","m7G","m7G","D","m7G"))
  mrs2 <- combineIntoModstrings(mrs, gr)
  expect_false(all(rs == c(mrs2[1:3],rs[4:5])))
  expect_true(all((rs == c(mrs2[1:3],rs[4:5]))[4:5]))
})
