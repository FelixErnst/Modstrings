#' @include Modstrings.R
#' @include Modstrings-ModString.R
#' @include Modstrings-ModStringSet.R
#' @include Modstrings-ModStringViews.R
NULL

# derived from Biostrings/R/matchPatter.R --------------------------------------
# .XString.matchPattern() and .XStringViews.matchPattern()
.ModString.matchPattern <- function(pattern,
                                    subject,
                                    max.mismatch,
                                    min.mismatch,
                                    with.indels,
                                    fixed,
                                    algorithm,
                                    count.only = FALSE){
  algo <- Biostrings:::normargAlgorithm(algorithm)
  if (Biostrings:::isCharacterAlgo(algo)){
    return(Biostrings:::.character.matchPattern(pattern,
                                                subject,
                                                max.mismatch,
                                                fixed,
                                                algo,
                                                count.only))
  }
  if (!is(subject, "ModString")){
    subject <- ModString(NULL, subject)
  }
  pattern <- .normargPattern(pattern, subject)
  max.mismatch <- Biostrings:::normargMaxMismatch(max.mismatch)
  min.mismatch <- Biostrings:::normargMinMismatch(min.mismatch, max.mismatch)
  with.indels <- Biostrings:::normargWithIndels(with.indels)
  fixed <- Biostrings:::normargFixed(fixed, subject)
  if (!assertive::is_a_bool(count.only)){
    stop("'count.only' must be TRUE or FALSE")
  }
  algo <- Biostrings:::selectAlgo(algo,
                                  pattern,
                                  max.mismatch,
                                  min.mismatch,
                                  with.indels,
                                  fixed)
  C_ans <- .Call2("XString_match_pattern",
                  pattern,
                  subject,
                  max.mismatch,
                  min.mismatch,
                  with.indels,
                  fixed,
                  algo,
                  count.only,
                  PACKAGE = "Biostrings")
  if(count.only){
    return(C_ans)
  }
  unsafe.newModStringViews(subject,
                           start(C_ans),
                           width(C_ans))
}

.ModStringViews.matchPattern <- function(pattern,
                                         subject,
                                         max.mismatch,
                                         min.mismatch,
                                         with.indels,
                                         fixed,
                                         algorithm,
                                         count.only = FALSE){
  algo <- Biostrings:::normargAlgorithm(algorithm)
  if (Biostrings:::isCharacterAlgo(algo)){
    stop("'subject' must be a single (non-empty) string ",
         "for this algorithm")
  }
  pattern <- .normargPattern(pattern, subject)
  max.mismatch <- Biostrings:::normargMaxMismatch(max.mismatch)
  min.mismatch <- Biostrings:::normargMinMismatch(min.mismatch, max.mismatch)
  with.indels <- Biostrings:::normargWithIndels(with.indels)
  fixed <- Biostrings:::normargFixed(fixed, subject)
  if (!assertive::is_a_bool(count.only)){
    stop("'count.only' must be TRUE or FALSE")
  }
  algo <- Biostrings:::selectAlgo(algo,
                                  pattern,
                                  max.mismatch,
                                  min.mismatch,
                                  with.indels,
                                  fixed)
  C_ans <- .Call2("XStringViews_match_pattern",
                  pattern,
                  subject(subject),
                  start(subject),
                  width(subject),
                  max.mismatch,
                  min.mismatch,
                  with.indels,
                  fixed,
                  algo,
                  count.only,
                  PACKAGE = "Biostrings")
  if (count.only){
    return(C_ans)
  }
  unsafe.newModStringViews(subject(subject),
                           start(C_ans),
                           width(C_ans))
}

setMethod("matchPattern", "ModString",
          function(pattern,
                   subject,
                   max.mismatch = 0,
                   min.mismatch = 0,
                   with.indels = FALSE,
                   fixed = TRUE,
                   algorithm = "auto")
            .ModString.matchPattern(pattern, subject, max.mismatch, min.mismatch,
                                    with.indels, fixed, algorithm)
)

setMethod("matchPattern", "ModStringViews",
          function(pattern,
                   subject,
                   max.mismatch = 0,
                   min.mismatch = 0,
                   with.indels = FALSE,
                   fixed = TRUE,
                   algorithm = "auto")
            .ModStringViews.matchPattern(pattern, subject, max.mismatch, 
                                         min.mismatch, with.indels, fixed,
                                         algorithm)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("matchPattern", "MaskedModString",
          function(pattern,
                   subject,
                   max.mismatch = 0,
                   min.mismatch = 0,
                   with.indels = FALSE,
                   fixed = TRUE,
                   algorithm = "auto")
            matchPattern(pattern,
                         toModStringViewsOrModString(subject),
                         max.mismatch = max.mismatch,
                         min.mismatch = min.mismatch,
                         with.indels = with.indels,
                         fixed = fixed,
                         algorithm = algorithm)
)


# derived from Biostrings/R/matchPattern.R -------------------------------------

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "ModString",
          function(pattern,
                   subject,
                   max.mismatch = 0,
                   min.mismatch = 0,
                   with.indels = FALSE,
                   fixed = TRUE,
                   algorithm = "auto")
            .ModString.matchPattern(pattern, subject, max.mismatch, min.mismatch,
                                    with.indels, fixed, algorithm, 
                                    count.only = TRUE)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "ModStringViews",
          function(pattern,
                   subject,
                   max.mismatch = 0,
                   min.mismatch = 0,
                   with.indels = FALSE,
                   fixed = TRUE,
                   algorithm = "auto")
            .ModStringViews.matchPattern(pattern, subject, max.mismatch, 
                                         min.mismatch, with.indels, fixed, 
                                         algorithm, count.only = TRUE)
)

### Dispatch on 'subject' (see signature of generic).
setMethod("countPattern", "MaskedModString",
          function(pattern,
                   subject,
                   max.mismatch = 0,
                   min.mismatch = 0,
                   with.indels = FALSE,
                   fixed = TRUE,
                   algorithm = "auto")
            countPattern(pattern,
                         toModStringViewsOrModString(subject),
                         max.mismatch = max.mismatch,
                         min.mismatch = min.mismatch,
                         with.indels = with.indels,
                         fixed = fixed,
                         algorithm = algorithm)
)

# derived from Biostrings/R/matchPattern.R -------------------------------------
# The "vmatchPattern" and "vcountPattern" methods.
#
# These are vectorized versions of matchPattern() and countPattern().
# vmatchPattern() returns an MIndex object and vcountPattern() an integer
# vector (like matchPDict() and countPDict(), respectively).
#

.ModStringSet.vmatchPattern <- function(pattern, 
                                        subject,
                                        max.mismatch,
                                        min.mismatch,
                                        with.indels,
                                        fixed,
                                        algorithm,
                                        count.only = FALSE)
{
  if (!assertive::is_a_bool(count.only)) {
    stop("'count.only' must be TRUE or FALSE")
  }
  if (!is(subject, "ModStringSet")){
    subject <- ModStringSet(NULL, subject)
  }
  algo <- Biostrings:::normargAlgorithm(algorithm)
  if (Biostrings:::isCharacterAlgo(algo)){
    stop("'subject' must be a single (non-empty) string ", 
         "for this algorithm")
  }
  # this converts the patter into ModString. Therefore we don't need to convert
  # anything into a one byte letter. this is done automaitcally. But we need
  # the function call to a modified .normargPattern()
  pattern <- .normargPattern(pattern, subject)
  max.mismatch <- Biostrings:::normargMaxMismatch(max.mismatch)
  min.mismatch <- Biostrings:::normargMinMismatch(min.mismatch, max.mismatch)
  with.indels <- Biostrings:::normargWithIndels(with.indels)
  fixed <- Biostrings:::normargFixed(fixed, subject)
  algo <- Biostrings:::selectAlgo(algo, 
                                  pattern,
                                  max.mismatch,
                                  min.mismatch,
                                  with.indels,
                                  fixed)
  # because MIndex objects do not support variable-width matches yet
  if (algo == "indels" && !count.only){
    stop("vmatchPattern() does not support indels yet")
  }
  C_ans <- .Call2("XStringSet_vmatch_pattern",
                  pattern,
                  subject,
                  max.mismatch,
                  min.mismatch,
                  with.indels,
                  fixed,
                  algo,
                  ifelse(count.only, "MATCHES_AS_COUNTS", "MATCHES_AS_ENDS"),
                  PACKAGE = "Biostrings")
  if (count.only){
    return(C_ans)
  }
  ans_width0 <- rep.int(length(pattern), length(subject))
  new("ByPos_MIndex",
      width0 = ans_width0,
      NAMES = names(subject),
      ends = C_ans)
}

setMethod("vmatchPattern","ModStringSet",
          function(pattern,
                   subject,
                   max.mismatch = 0,
                   min.mismatch = 0,
                   with.indels = FALSE,
                   fixed = TRUE,
                   algorithm = "auto")
            .ModStringSet.vmatchPattern(pattern, subject, max.mismatch,
                                        min.mismatch, with.indels, fixed,
                                        algorithm)
)

setMethod("vcountPattern", "ModStringSet",
          function(pattern,
                   subject,
                   max.mismatch = 0,
                   min.mismatch = 0,
                   with.indels = FALSE,
                   fixed = TRUE,
                   algorithm = "auto")
            .ModStringSet.vmatchPattern(pattern, subject, max.mismatch,
                                        min.mismatch, with.indels, fixed,
                                        algorithm, count.only = TRUE)
)

setMethod("vcountPattern", "ModStringViews",
          function(pattern,
                   subject,
                   max.mismatch = 0,
                   min.mismatch = 0,
                   with.indels = FALSE,
                   fixed = TRUE,
                   algorithm = "auto")
            vcountPattern(pattern,
                          toModStringViewsOrModString(subject),
                          max.mismatch = max.mismatch,
                          min.mismatch = min.mismatch,
                          with.indels = with.indels,
                          fixed = fixed,
                          algorithm = algorithm)
)
