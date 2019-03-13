#' @include Modstrings.R
NULL

#' @name letterFrequency
#' @aliases hasOnlyBaseLetters hasOnlyBaseLetters consensusString 
#' alphabetFrequency letterFrequencyInSlidingView
#' hasOnlyBaseLetters,ModDNAString-method
#' hasOnlyBaseLetters,ModDNAStringSet-method
#' hasOnlyBaseLetters,ModRNAString-method
#' hasOnlyBaseLetters,ModRNAStringSet-method
#' 
#' @title Calculate the frequency of letters in nucleotide sequence with 
#' modifications, or the consensus matrix of a set of sequences
#' 
#' @description 
#' These functions follow the same principle as the 
#' \code{\link[Biostrings:letterFrequency]{Biostrings}} functions. Please be
#' aware, that the matices can become quite large, since the alphabet of 
#' \code{ModString} objects contains more letters.
#'
#' @param x a \code{\link{ModString}}, a \code{\link{ModStringSet}}, 
#' a \code{\link{ModStringViews}} or a \code{\link{MaskedModString}} object.
#' @param as.prob \code{TRUE} or \code{FALSE} (default): Should the result be 
#' returned as probabilities instead of counts? (\code{sum per column = 1})
#' @param baseOnly \code{TRUE} or \code{FALSE} (default): Should the result omit
#' occurances of the letters \code{N.-+}?
#' @param letters See \code{\link[Biostrings]{letterFrequency}}.
#' @param OR See \code{\link[Biostrings]{letterFrequency}}.
#' @param threshold Since the amiguityMap is fixed to \code{"?"} for 
#' \code{ModString} objects, only the treshold can be set (default 
#' \code{threshold = 0.25})
#' @param collapse \code{TRUE} or \code{FALSE} (default): Should the results 
#' summed up all elements for \code{ModStringSet} or \code{ModStringViews} 
#' objects or reported per element.
#' @param shift See \code{\link[Biostrings]{letterFrequency}}.
#' @param width See \code{\link[Biostrings]{letterFrequency}}.
#' @param ... See \code{\link[Biostrings]{letterFrequency}}.
#' 
#' @return a matrix with the results (letter x pos).
#' 
#' @export 
#'
#' @examples
#' mod <- ModDNAString(paste(alphabet(ModDNAString()), collapse = ""))
#' hasOnlyBaseLetters(mod)
#' alphabetFrequency(mod)
NULL

#' @rdname letterFrequency
#' @export
setMethod("hasOnlyBaseLetters", "ModDNAString",
          function(x) hasOnlyBaseLetters(DNAString(x))
)
#' @rdname letterFrequency
#' @export
setMethod("hasOnlyBaseLetters", "ModRNAString",
          function(x) hasOnlyBaseLetters(RNAString(x))
)
#' @export
setMethod("hasOnlyBaseLetters", "ModDNAStringSet",
          function(x) hasOnlyBaseLetters(DNAStringSet(x))
)
#' @export
setMethod("hasOnlyBaseLetters", "ModRNAStringSet",
          function(x) hasOnlyBaseLetters(RNAStringSet(x))
)

#' @rdname letterFrequency
#' @export
setMethod("alphabetFrequency", "ModDNAString",
          function(x, as.prob = FALSE, baseOnly = FALSE)
          {
            ans <- .XString.nucleotide_frequency(x, as.prob, baseOnly)
            names(ans) <- names(xscodes(x, baseOnly = baseOnly,
                                        multiByteLetterNames = TRUE))
            ans
          }
)
#' @rdname letterFrequency
#' @export
setMethod("alphabetFrequency", "ModRNAString",
          function(x, as.prob = FALSE, baseOnly = FALSE)
          {
            ans <- .XString.nucleotide_frequency(x, as.prob, baseOnly)
            names(ans) <- names(xscodes(x, baseOnly = baseOnly,
                                        multiByteLetterNames = TRUE))
            ans
          }
)
#' @rdname letterFrequency
#' @export
setMethod("alphabetFrequency", "ModDNAStringSet",
          function(x, as.prob = FALSE, collapse = FALSE, baseOnly = FALSE)
          {
            ans <- .XStringSet.nucleotide_frequency(x, as.prob, collapse,
                                                    baseOnly)
            colnames(ans) <- names(xscodes(x, baseOnly = baseOnly,
                                           multiByteLetterNames = TRUE))
            ans
          }
)
#' @rdname letterFrequency
#' @export
setMethod("alphabetFrequency", "ModRNAStringSet",
          function(x, as.prob = FALSE, collapse = FALSE, baseOnly = FALSE)
          {
            ans <- .XStringSet.nucleotide_frequency(x, as.prob, collapse,
                                                    baseOnly)
            colnames(ans) <- names(xscodes(x, baseOnly = baseOnly,
                                           multiByteLetterNames = TRUE))
            ans
          }
)

#' @rdname letterFrequency
#' @export
setMethod("alphabetFrequency", "MaskedModString",
          function(x, as.prob = FALSE, ...){
            y <- as(x, "ModStringViews")
            alphabetFrequency(y, as.prob = as.prob, collapse = TRUE, ...)
          }
)

#' @rdname letterFrequency
#' @export
setMethod("letterFrequency", "ModStringViews",
          function(x, letters, OR = "|", as.prob = FALSE,  ...)
            letterFrequency(as(x, "ModStringSet"), letters = letters, OR = OR,
                            as.prob = as.prob, ...)
)
#' @rdname letterFrequency
#' @export
setMethod("letterFrequency", "MaskedModString",
          function(x, letters, OR = "|", as.prob = FALSE)
            letterFrequency(as(x, "ModStringViews"), letters = letters, OR = OR,
                            as.prob = as.prob, collapse = TRUE)
)


# derived from Biostrings/R/letterFrequency.R ----------------------------------
#' @rdname letterFrequency
#' @export
setMethod(
  "consensusMatrix", "ModStringSet",
  function(x, as.prob = FALSE, shift = 0L, width = NULL, baseOnly = FALSE)
  {
    ans <- callNextMethod()
    rownames(ans) <- names(xscodes(x, multiByteLetterNames = TRUE))
    ans
  }
)

#' @rdname letterFrequency
#' @export
setMethod("consensusString", "ModDNAStringSet",
          function(x, threshold = 0.25, shift = 0L, width = NULL)
          {
            consensusString(consensusMatrix(x, as.prob = TRUE, shift = shift,
                                            width = width),
                            ambiguityMap = "?",
                            threshold = threshold)
          }
)

#' @rdname letterFrequency
#' @export
setMethod("consensusString", "ModRNAStringSet",
          function(x, threshold = 0.25, shift = 0L, width = NULL)
          {
            consensusString(consensusMatrix(x, as.prob = TRUE, shift = shift,
                                            width = width),
                            ambiguityMap = "?",
                            threshold = threshold)
          }
)

#' @rdname letterFrequency
#' @export
setMethod("consensusString", "ModStringViews",
          function(x, threshold, shift = 0L, width = NULL)
          {
            x <- as(x, "ModStringSet")
            callGeneric()
          }
)
