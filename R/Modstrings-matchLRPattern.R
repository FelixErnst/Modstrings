#' @include Modstrings.R
#' @include Modstrings-ModString.R
#' @include Modstrings-ModStringSet.R
#' @include Modstrings-ModStringViews.R
#' @include Modstrings-matchPattern.R
NULL

# derived from Biostrings/R/matchLRPattern.R -----------------------------------

setMethod("matchLRPatterns", "ModString", 
          function(Lpattern,
                   Rpattern,
                   max.gaplength,
                   subject,
                   max.Lmismatch = 0,
                   max.Rmismatch = 0,
                   with.Lindels = FALSE,
                   with.Rindels = FALSE,
                   Lfixed = TRUE,
                   Rfixed = TRUE){
            ans_start <- ans_end <- integer(0)
            Lmatches <- matchPattern(Lpattern,
                                     subject,
                                     max.mismatch = max.Lmismatch,
                                     with.indels = with.Lindels,
                                     fixed = Lfixed)
            if (length(Lmatches) != 0L) {
              Rmatches <- matchPattern(Rpattern,
                                       subject,
                                       max.mismatch = max.Rmismatch,
                                       with.indels = with.Rindels,
                                       fixed = Rfixed)
              if (length(Rmatches) != 0L) {
                for (i in seq_len(length(Lmatches))) {
                  gaplength <- start(Rmatches) - end(Lmatches)[i] - 1L
                  jj <- which(0L <= gaplength & gaplength <= max.gaplength)
                  ans_start <- c(ans_start,
                                 rep.int(start(Lmatches)[i],
                                         length(jj)))
                  ans_end <- c(ans_end, end(Rmatches)[jj])
                }
              }
            }
            ans_width <- ans_end - ans_start + 1L
            unsafe.newModStringViews(subject, ans_start, ans_width)
          }
)

setMethod("matchLRPatterns", "ModStringViews",
          function(Lpattern,
                   Rpattern,
                   max.gaplength,
                   subject,
                   max.Lmismatch = 0,
                   max.Rmismatch = 0,
                   with.Lindels = FALSE,
                   with.Rindels = FALSE,
                   Lfixed = TRUE,
                   Rfixed = TRUE){
            ans_start <- ans_width <- integer(0)
            subject <- trim(subject)
            Lcounts <- vcountPattern(Lpattern,
                                     subject,
                                     max.mismatch = max.Lmismatch,
                                     with.indels = with.Lindels,
                                     fixed = Lfixed)
            candidates <- which(Lcounts > 0L)
            if (length(candidates) != 0L) {
              subject <- subject[candidates]
              Lcounts <- Lcounts[candidates]
              Lmatches <- matchPattern(Lpattern,
                                       subject,
                                       max.mismatch = max.Lmismatch,
                                       with.indels = with.Lindels,
                                       fixed = Lfixed)
              Rranges <- IRanges::IRanges(start = end(Lmatches) + 1L,
                                          width = (max.gaplength + 
                                            nchar(Rpattern)))
              Rranges <- pintersect(Rranges,
                                    rep(as(subject, "IRanges"), Lcounts))
              Rsubject <- as(Views(subject(subject), Rranges), "XStringSet")
              Rmatches <- vmatchPattern(Rpattern,
                                        Rsubject,
                                        max.mismatch = max.Rmismatch,
                                        with.indels = with.Rindels,
                                        fixed = Rfixed)
              Rcounts <- elementNROWS(Rmatches)
              Roffset <- unlist(endIndex(Rmatches), use.names=FALSE)
              if (length(Roffset) != 0L) {
                ans_start <- rep.int(start(Lmatches), Rcounts)
                ans_width <- rep.int(width(Lmatches), Rcounts) + Roffset
              }
            }
            unsafe.newModStringViews(subject(subject),
                                     ans_start,
                                     ans_width)
          }
)

setMethod("matchLRPatterns", "MaskedModString",
          function(Lpattern,
                   Rpattern,
                   max.gaplength,
                   subject,
                   max.Lmismatch = 0,
                   max.Rmismatch = 0,
                   with.Lindels = FALSE,
                   with.Rindels = FALSE,
                   Lfixed = TRUE,
                   Rfixed = TRUE)
            matchLRPatterns(Lpattern,
                            Rpattern,
                            max.gaplength,
                            toModStringViewsOrModString(subject),
                            max.Lmismatch = max.Lmismatch,
                            max.Rmismatch = max.Rmismatch,
                            with.Lindels = with.Lindels,
                            with.Rindels = with.Rindels,
                            Lfixed = Lfixed,
                            Rfixed = Rfixed)
)

