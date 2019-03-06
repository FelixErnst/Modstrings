#' @include Modstrings.R
NULL

# Modstrings-ModString.R -------------------------------------------------------
.call_new_CHARACTER_from_XString <- function(x){
  .Call2("new_CHARACTER_from_XString",
         x, 
         NULL,
         PACKAGE = "Biostrings")
}
.call_new_XString_from_CHARACTER <- function(classname, x, start, width){
  .Call2("new_XString_from_CHARACTER",
         classname,
         x,
         start,
         width,
         NULL,
         PACKAGE = "Biostrings")
}

# Modstrings-ModStringSet.R ----------------------------------------------------
.call_new_CHARACTER_from_XStringSet <- function(x){
  .Call2("new_CHARACTER_from_XStringSet",
         x, 
         NULL,
         PACKAGE = "Biostrings")
}
.call_new_XStringSet_from_CHARACTER <- function(classname, elementType, x, 
                                                start, width){
  .Call2("new_XStringSet_from_CHARACTER",
         classname,
         elementType,
         x,
         start,
         width,
         NULL,
         PACKAGE = "Biostrings")
}

# Modstrings-letterFrequency.R -------------------------------------------------
# 
.call_XString_letter_frequency <- function(x, codes, baseOnly){
  .Call2("XString_letter_frequency",
         x,
         codes,
         baseOnly,
         PACKAGE = "Biostrings")
}

# Don not confuse with the next .call function. Very similiar names
.call_XStringSet_letter_frequency <- function(x, collapse, codes, baseOnly){
  .Call2("XStringSet_letter_frequency",
         x, 
         collapse,
         codes,
         baseOnly,
         PACKAGE = "Biostrings")
}
# See comment above
.call_XStringSet_letterFrequency <- function(x, single_codes, colmap, colnames,
                                               collapse){
  .Call2("XStringSet_letterFrequency",
         x,
         single_codes,
         colmap,
         colnames,
         collapse,
         PACKAGE = "Biostrings")
}

.call_XString_letterFrequencyInSlidingView <- function(x, view.width, 
                                                       single_codes, colmap,
                                                       colnames){
  .Call2("XString_letterFrequencyInSlidingView",
         x,
         view.width,
         single_codes,
         colmap,
         colnames,
         PACKAGE="Biostrings")
}

.call_XStringSet_consensus_matrix <- function(x, shift, width, baseOnly,
                                              codes){
  .Call2("XStringSet_consensus_matrix",
         x,
         shift,
         width,
         baseOnly,
         codes,
         PACKAGE="Biostrings")
}

# Modstrings-matching-lowlevel.R -----------------------------------------------
.call_XString_match_pattern_at <- function(pattern, subject, at, at.type,
                                           max.mismatch, min.mismatch,
                                           with.indels, fixed, ans.type,
                                           auto.reduce.pattern){
  .Call2("XString_match_pattern_at",
         pattern,
         subject,
         at,
         at.type,
         max.mismatch,
         min.mismatch,
         with.indels,
         fixed,
         ans.type,
         auto.reduce.pattern,
         PACKAGE="Biostrings")
}
.call_XStringSet_vmatch_pattern_at <- function(pattern, subject, at, at.type,
                                               max.mismatch, min.mismatch,
                                               with.indels, fixed, ans.type,
                                               auto.reduce.pattern){
  .Call2("XStringSet_vmatch_pattern_at",
         pattern,
         subject,
         at,
         at.type,
         max.mismatch,
         min.mismatch,
         with.indels,
         fixed,
         ans.type,
         auto.reduce.pattern,
         PACKAGE="Biostrings")
}

# Modstrings-matchPattern.R ----------------------------------------------------
.call_XString_match_pattern <- function(pattern, subject, max.mismatch,
                                        min.mismatch, with.indels, fixed,
                                        algo, count.only){
  .Call2("XString_match_pattern",
         pattern,
         subject,
         max.mismatch,
         min.mismatch,
         with.indels,
         fixed,
         algo,
         count.only,
         PACKAGE = "Biostrings")
}

.call_XStringViews_match_pattern <- function(pattern, subject, start, width,
                                             max.mismatch, min.mismatch,
                                             with.indels, fixed, algo, 
                                             count.only){
  .Call2("XStringViews_match_pattern",
         pattern,
         subject,
         start,
         width,
         max.mismatch,
         min.mismatch,
         with.indels,
         fixed,
         algo,
         count.only,
         PACKAGE = "Biostrings")
}

.call_XStringSet_vmatch_pattern <- function(pattern, subject, max.mismatch,
                                            min.mismatch, with.indels, fixed,
                                            algo, count_only){
  .Call2("XStringSet_vmatch_pattern",
         pattern,
         subject,
         max.mismatch,
         min.mismatch,
         with.indels,
         fixed,
         algo,
         count_only,
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
