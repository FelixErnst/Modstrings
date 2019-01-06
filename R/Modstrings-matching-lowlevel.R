#' @include Modstrings.R
NULL

# derived from Biostrings/R/lowlevel-matching.R --------------------------------
# Helper functions (not exported) used by matching functions from other
# files (like matchPattern(), matchPDict(), etc...) to check and normalize
# their arguments.

.normargSubject <- function(subject,
                            argname = "subject"){
  if (is(subject, "ModString") || is(subject, "ModStringSet")){
    return(subject)
  }
  if (!is.character(subject)){
    stop("'", argname, "' must be a character vector, ",
         "or an ModString or ModStringSet object")
  }
  if (length(subject) == 1L){
    subject <- as(subject, "ModString")
  } else {
    subject <- as(subject, "ModStringSet")
  }
  subject
}

.normargPattern <- function(pattern,
                            subject,
                            argname = "pattern"){
  subject_baseclass <- Biostrings:::xsbaseclass(subject)
  if (is(pattern, "ModString")) {
    if (base_class_name(pattern) == subject_baseclass){
      return(pattern)
    }
  } else if (!assertive::is_a_non_empty_string(pattern)){
    stop("'", argname, "' must be a single string or an ModString object")
  }
  pattern <- try(ModString(seqtype(subject), pattern))
  if (is(pattern, "try-error")){
    stop("could not turn '", argname, "' into a ",
         subject_baseclass, " instance")
  }
  pattern
}

# derived from Biostrings/R/lowlevel-matching.R --------------------------------
# needed because of calls to
# - .normargSubject
# - .normargPattern


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# .matchPatternAt()
#
# This is the horse-power behind the neditStartingAt(), neditEndingAt(),
# isMatchingStartingAt(), isMatchingEndingAt(), which.isMatchingStartingAt()
# and which.isMatchingEndingAt() low-level matching functions defined later
# in this file.
# If 'at.type == 0' then 'at' contains starting positions, otherwise it
# contains ending positions.

.matchPatternAt <- function(pattern,
                            subject,
                            at,
                            at.type,
                            max.mismatch,
                            min.mismatch,
                            with.indels,
                            fixed,
                            ans.type,
                            auto.reduce.pattern = FALSE){
  subject <- .normargSubject(subject)
  pattern <- .normargPattern(pattern, subject)
  if (!is.numeric(at)) {
    what <- if (at.type == 0L) "starting.at" else "ending.at"
    stop("'", what, "' must be a vector of integers")
  }
  if (!is.integer(at))
    at <- as.integer(at)
  
  if (auto.reduce.pattern) {
    at.length <- length(at)
    P <- nchar(pattern)
    if (at.length == 1)
      at <- rep.int(at, P)
    else if (at.length != P || length(unique(at)) > 1)
      stop("With 'auto.reduce.pattern', 'at' must be a single integer")
  }
  
  if (ans.type == 0L) {
    max.mismatch <- length(pattern)
  } else {
    if (!is.numeric(max.mismatch))
      stop("'max.mismatch' must be a vector of integers")
    if (!is.integer(max.mismatch))
      max.mismatch <- as.integer(max.mismatch)
    if (!is.numeric(min.mismatch))
      stop("'min.mismatch' must be a vector of integers")
    if (!is.integer(min.mismatch))
      min.mismatch <- as.integer(min.mismatch)
  }
  with.indels <- Biostrings:::normargWithIndels(with.indels)
  fixed <- Biostrings:::normargFixed(fixed, subject)
  if (is(subject, "ModString")){
    ans <- .Call2("XString_match_pattern_at",
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
  } else {
    ans <- .Call2("XStringSet_vmatch_pattern_at",
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
  ans
}


# derived from Biostrings/R/lowlevel-matching.R --------------------------------

setMethod("neditStartingAt", "ModString",
          function(pattern,
                   subject,
                   starting.at = 1,
                   with.indels = FALSE,
                   fixed = TRUE)
            .matchPatternAt(pattern,
                            subject,
                            starting.at,
                            0L,
                            NA,
                            NA,
                            with.indels,
                            fixed,
                            0L)
)

setMethod("neditStartingAt", "ModStringSet",
          function(pattern,
                   subject,
                   starting.at = 1,
                   with.indels = FALSE,
                   fixed = TRUE)
            .matchPatternAt(pattern,
                            subject,
                            starting.at,
                            0L,
                            NA,
                            NA,
                            with.indels,
                            fixed,
                            0L)
)

setMethod("neditEndingAt", "ModString",
          function(pattern,
                   subject,
                   ending.at = 1,
                   with.indels = FALSE,
                   fixed = TRUE)
            .matchPatternAt(pattern,
                            subject,
                            ending.at,
                            1L,
                            NA,
                            NA,
                            with.indels,
                            fixed,
                            0L)
)

setMethod("neditEndingAt", "ModStringSet",
          function(pattern,
                   subject,
                   ending.at = 1,
                   with.indels = FALSE,
                   fixed = TRUE)
            .matchPatternAt(pattern,
                            subject,
                            ending.at,
                            1L,
                            NA,
                            NA,
                            with.indels,
                            fixed,
                            0L)
)


# derived from Biostrings/R/lowlevel-matching.R --------------------------------

# isMatchingStartingAt(), isMatchingEndingAt() and isMatchingAt().
#
# 'starting.at' (or 'ending.at') must be integer vectors containing the
# starting (or ending) positions of the pattern relatively to the subject.
# These functions return a logical vector of the same length as
# 'starting.at' (or 'ending.at').

setMethod("isMatchingStartingAt", "ModString",
          function(pattern,
                   subject,
                   starting.at = 1,
                   max.mismatch = 0,
                   min.mismatch = 0,
                   with.indels = FALSE,
                   fixed = TRUE)
            .matchPatternAt(pattern,
                            subject,
                            starting.at,
                            0L,
                            max.mismatch,
                            min.mismatch,
                            with.indels,
                            fixed,
                            1L)
)

setMethod("isMatchingStartingAt", "ModStringSet",
          function(pattern,
                   subject,
                   starting.at = 1,
                   max.mismatch = 0,
                   min.mismatch = 0,
                   with.indels = FALSE,
                   fixed = TRUE)
            .matchPatternAt(pattern,
                            subject,
                            starting.at,
                            0L,
                            max.mismatch,
                            min.mismatch,
                            with.indels,
                            fixed,
                            1L)
)

setMethod("isMatchingEndingAt", "ModString",
          function(pattern,
                   subject,
                   ending.at = 1,
                   max.mismatch = 0,
                   min.mismatch = 0,
                   with.indels = FALSE,
                   fixed = TRUE)
            .matchPatternAt(pattern,
                            subject,
                            ending.at,
                            0L,
                            max.mismatch,
                            min.mismatch,
                            with.indels,
                            fixed,
                            1L)
)

setMethod("isMatchingEndingAt", "ModStringSet",
          function(pattern,
                   subject,
                   ending.at = 1,
                   max.mismatch = 0,
                   min.mismatch = 0,
                   with.indels = FALSE,
                   fixed = TRUE)
            .matchPatternAt(pattern,
                            subject,
                            ending.at,
                            0L,
                            max.mismatch,
                            min.mismatch,
                            with.indels,
                            fixed,
                            1L)
)


# derived from Biostrings/R/lowlevel-matching.R --------------------------------
# which.isMatchingStartingAt(), which.isMatchingEndingAt() and
# which.isMatchingAt().
#
# 'starting.at' (or 'ending.at') must be integer vectors containing the
# starting (or ending) positions of the pattern relatively to the subject.
# These functions return the lowest *index* in 'starting.at' (or 'ending.at')
# for which a match occurred (or NA if no match occurred).


setMethod("which.isMatchingStartingAt", "ModString",
          function(pattern,
                   subject,
                   starting.at = 1,
                   max.mismatch = 0,
                   min.mismatch = 0,
                   with.indels = FALSE,
                   fixed = TRUE,
                   follow.index = FALSE,
                   auto.reduce.pattern = FALSE)
            .matchPatternAt(pattern,
                            subject,
                            starting.at,
                            0L,
                            max.mismatch,
                            min.mismatch,
                            with.indels,
                            fixed,
                            Biostrings:::.to.ans.type(follow.index),
                            auto.reduce.pattern)
)

setMethod("which.isMatchingStartingAt", "ModStringSet",
          function(pattern,
                   subject,
                   starting.at = 1,
                   max.mismatch = 0,
                   min.mismatch = 0,
                   with.indels = FALSE,
                   fixed = TRUE,
                   follow.index = FALSE,
                   auto.reduce.pattern = FALSE)
            .matchPatternAt(pattern,
                            subject,
                            starting.at,
                            0L,
                            max.mismatch,
                            min.mismatch,
                            with.indels,
                            fixed,
                            Biostrings:::.to.ans.type(follow.index),
                            auto.reduce.pattern)
)

setMethod("which.isMatchingEndingAt", "ModString",
          function(pattern,
                   subject,
                   ending.at = 1,
                   max.mismatch = 0,
                   min.mismatch = 0,
                   with.indels = FALSE,
                   fixed = TRUE,
                   follow.index = FALSE,
                   auto.reduce.pattern = FALSE)
            .matchPatternAt(pattern,
                            subject,
                            ending.at,
                            1L,
                            max.mismatch,
                            min.mismatch,
                            with.indels,
                            fixed,
                            Biostrings:::.to.ans.type(follow.index),
                            auto.reduce.pattern)
)

setMethod("which.isMatchingEndingAt", "ModStringSet",
          function(pattern,
                   subject,
                   ending.at = 1,
                   max.mismatch = 0,
                   min.mismatch = 0,
                   with.indels = FALSE,
                   fixed = TRUE,
                   follow.index = FALSE,
                   auto.reduce.pattern = FALSE)
            .matchPatternAt(pattern,
                            subject,
                            ending.at,
                            1L,
                            max.mismatch,
                            min.mismatch,
                            with.indels,
                            fixed,
                            Biostrings:::.to.ans.type(follow.index),
                            auto.reduce.pattern)
)


# derived from Biostrings/R/lowlevel-matching.R --------------------------------
# hasLetterAt()

hasLetterAt <- function(x, letter, at, fixed=TRUE){
  if (!is(x, "ModStringSet")) {
    if (!is.character(x) && !is(x, "ModString")){
      stop("'x' must be a character vector, or an ModString or ModStringSet ",
           "object")
    }
    x <- as(x, "ModStringSet")
  }
  if (!is(letter, "ModString")) {
    if (!assertive::is_a_non_empty_string(letter)){
      stop("'letter' must be a character string or an ModString object")
    }
    letter <- ModString(seqtype(x), letter)
  } else {
    if (seqtype(letter) != seqtype(x)){
      stop("'x' and 'letter' must have the same sequence type")
    }
  }
  if (!is.numeric(at)){
    stop("'at' must be a vector of integers")
  }
  if (length(at) != length(letter)){
    stop("'letter' and 'at' must have the same length")
  }
  if (!is.integer(at)){
    at <- as.integer(at)
  }
  if (any(is.na(at))){
    stop("'at' cannot contain NAs")
  }
  fixed <- Biostrings:::normargFixed(fixed, x)
  
  .hasLetterAt1 <- function(x, l1, at1){
    ans <- .Call2("XStringSet_vmatch_pattern_at",
                  l1,
                  x,
                  at1,
                  0L,
                  0L,
                  0L,
                  FALSE,
                  fixed,
                  1L,
                  FALSE,
                  PACKAGE="Biostrings")
    ans[at1 < 1 | at1 > width(x)] <- NA
    ans
  }
  vapply(seq_len(length(letter)),
         function(i)
           .hasLetterAt1(x, subseq(letter, start = i, width = 1L), at[i]),
         logical(1))
}

