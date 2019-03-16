#' @include Modstrings.R
NULL

# Modstrings-ModStringSet.R ----------------------------------------------------
.call_new_CHARACTER_from_XStringSet <- function(x){
  .Call2("new_CHARACTER_from_XStringSet",
         x, 
         NULL,
         PACKAGE = "Biostrings")
}

# Modstrings-modifyNucleotide.R ------------------------------------------------
.call_XString_replace_letter_at <- function(x, at, letter, verbose){
  .Call2("XString_replace_letter_at",
         x,
         at,
         letter,
         NULL,
         if.not.extending = "replace",
         verbose,
         PACKAGE = "Biostrings")
}

# Modstrings-ModStringSet-io.R -------------------------------------------------
.call_fasta_index <- function(filexp_list,
                              nrec,
                              skip,
                              seek.first.rec){
  # Guess what!? it does not like invalid one-letter sequence codes
  # suppress the warning. Results are valid noneheless for this purpose
  # However the offset cannot be used, since readBin output cannot be encoded
  # correctly
  suppressWarnings(C_ans <- .Call2("fasta_index",
                                   filexp_list,
                                   nrec,
                                   skip,
                                   seek.first.rec,
                                   NULL,
                                   PACKAGE = "Biostrings"))
  C_ans
}
