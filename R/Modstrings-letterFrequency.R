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
#' \code{ModString} obejcts contains more letters.
#'
#' @param x a \code{\link{ModString}} or a \code{\link{ModStringSet}} object.
#' @param as.prob \code{TRUE} or \code{FALSE} (default): Should the result be 
#' returned as probabilities instead of counts? (\code{sum per column = 1})
#' @param baseOnly \code{TRUE} or \code{FALSE} (default): Should the result omit
#' occurances of the letters \code{N.-+}?
#' @param view.width For \code{letterFrequencyInSlidingView}, the constant 
#' (e.g. 35, 48, 1000) size of the "window" to slide along x. The specified 
#' letters are tabulated in each window of length view.width. The rows of the 
#' result (see value) correspond to the various windows.
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

# derived from Biostrings/R/letterFrequency.R ----------------------------------

.ModString.nucleotide_frequency <- function(x,
                                            as.prob,
                                            baseOnly){
  if (!assertive::is_a_bool(as.prob)){
    stop("'as.prob' must be TRUE or FALSE", call. = FALSE)
  }
  codes <- modscodes(seqtype(x),
                     baseOnly = baseOnly)
  codes[] <- as.integer(unlist(lapply(names(codes),charToRaw)))
  ans <- .Call2("XString_letter_frequency",
                x,
                codes,
                baseOnly,
                PACKAGE = "Biostrings")
  names(ans) <- names(modscodes(seqtype(x),
                                lettersOnly = TRUE,
                                baseOnly = baseOnly))
  if (as.prob){
    ans <- ans / nchar(x) # nchar(x) is sum(ans) but faster
  }
  ans
}

.ModStringSet.nucleotide_frequency <- function(x,
                                               as.prob,
                                               collapse,
                                               baseOnly){
  if (!assertive::is_a_bool(as.prob)){
    stop("'as.prob' must be TRUE or FALSE", call. = FALSE)
  }
  collapse <- Biostrings:::.normargCollapse(collapse)
  codes <- modscodes(seqtype(x),
                     baseOnly = baseOnly)
  codes[] <- as.integer(unlist(lapply(names(codes),charToRaw)))
  ans <- .Call2("XStringSet_letter_frequency",
                x, 
                collapse,
                codes,
                baseOnly,
                PACKAGE = "Biostrings")
  if (collapse) {
    names(ans) <- names(modscodes(seqtype(x),
                                  lettersOnly = TRUE,
                                  baseOnly = baseOnly))
    fun <- "sum"
  } else {
    colnames(ans) <- names(modscodes(seqtype(x),
                                     lettersOnly = TRUE,
                                     baseOnly = baseOnly))
    fun <- "nchar"
  }
  if (as.prob) {
    ans <- ans / do.call(fun,list(ans))
  }
  ans
}

#' @rdname letterFrequency
#' @export
setMethod("alphabetFrequency", "ModDNAString",
          function(x, as.prob = FALSE, baseOnly = FALSE)
            .ModString.nucleotide_frequency(x, as.prob, baseOnly)
)
#' @rdname letterFrequency
#' @export
setMethod("alphabetFrequency", "ModRNAString",
          function(x, as.prob = FALSE, baseOnly = FALSE)
            .ModString.nucleotide_frequency(x, as.prob, baseOnly)
)
#' @rdname letterFrequency
#' @export
setMethod("alphabetFrequency", "ModDNAStringSet",
          function(x, as.prob = FALSE, collapse = FALSE, baseOnly = FALSE)
            .ModStringSet.nucleotide_frequency(x, as.prob, collapse, baseOnly)
)
#' @rdname letterFrequency
#' @export
setMethod("alphabetFrequency", "ModRNAStringSet",
          function(x, as.prob = FALSE, collapse = FALSE, baseOnly = FALSE)
            .ModStringSet.nucleotide_frequency(x, as.prob, collapse, baseOnly)
)
#' @rdname letterFrequency
#' @export
setMethod("alphabetFrequency", "ModStringViews",
          function(x, as.prob = FALSE, ...){
            y <- Biostrings:::fromXStringViewsToStringSet(x)
            alphabetFrequency(y, as.prob = as.prob, ...)
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


# derived from Biostrings/R/letterFrequency.R ----------------------------------

.letterFrequency <- function(x,
                             view.width,
                             letters,
                             OR,
                             collapse = FALSE){
  ## letterFrequency / letterFrequencyInSlidingView switch
  is_sliding <- !is.na(view.width)
  single_letters <- Biostrings:::.normargLetters(letters, alphabet(x))
  OR <- Biostrings:::.normargOR(OR)
  codes <- modscodes(seqtype(x))
  codes[] <- as.integer(unlist(lapply(names(codes),charToRaw)))
  names(codes) <- names(modscodes(seqtype(x),
                                  lettersOnly = TRUE))
  single_codes <- codes[single_letters]
  ## Unless 'OR == 0', letters in multi-character elements of
  ## 'letters' are to be grouped (i.e. tabulated in common).
  ## We send a vector indicating the column (1-based) into which each
  ## letter in 'letters' should be tabulated.  For example, for
  ## 'letters = c("CG", "AT")' and 'OR != 0', we send 'c(1,1,2,2)'.
  ## The columns of the result are named accordingly using the OR symbol.
  nc <- nchar(letters)
  if (all(nc == 1L) || OR == 0) {
    colmap <- NULL
    colnames <- single_letters
  } else {
    colmap <- rep.int(seq_len(length(letters)), nc)
    colnames <- vapply(strsplit(letters, NULL, fixed=TRUE),
                       function(z){
                         paste(z, collapse = OR)
                       },
                       character(1))
  }
  if (is_sliding){
    ans <- .Call2("XString_letterFrequencyInSlidingView",
                  x,
                  view.width,
                  single_codes,
                  colmap,
                  colnames,
                  PACKAGE="Biostrings")
  } else {
    ans <- .Call2("XStringSet_letterFrequency",
                  x,
                  single_codes,
                  colmap,
                  colnames,
                  collapse,
                  PACKAGE="Biostrings")
  }
  ans
}

### Ensure 'view.width' is not NA
#' @rdname letterFrequency
#' @export
setMethod("letterFrequencyInSlidingView", "ModString",
          function(x,
                   view.width,
                   letters,
                   OR = "|",
                   as.prob = FALSE){
            view.width <- Biostrings:::.normargWidth(view.width,
                                                     "view.width")
            if (!assertive::is_a_bool(as.prob)){
              stop("'as.prob' must be TRUE or FALSE")
            }
            ans <- .letterFrequency(x,
                                    view.width,
                                    letters = letters,
                                    OR = OR)
            if (as.prob){
              ans <- ans / view.width
            }
            ans
          }
)
#' @rdname letterFrequency
#' @export
setMethod("letterFrequency", "ModStringSet",
          function(x,
                   letters,
                   OR="|",
                   as.prob = FALSE,
                   collapse = FALSE){
            if (!assertive::is_a_bool(as.prob)){
              stop("'as.prob' must be TRUE or FALSE")
            }
            if (!assertive::is_a_bool(collapse)){
              stop("'collapse' must be TRUE or FALSE")
            }
            ans <- .letterFrequency(x,
                                    NA,
                                    letters = letters,
                                    OR = OR,
                                    collapse = collapse)
            if (as.prob) {
              nc <- nchar(x)
              if (collapse){
                nc <- sum(nc)
              }
              ans <- ans / nc
            }
            ans
          }
)
#' @rdname letterFrequency
#' @export
setMethod("letterFrequency", "ModStringViews",
          function(x,
                   letters,
                   OR = "|",
                   as.prob = FALSE,
                   ...)
            letterFrequency(as(x,
                               "ModStringSet"),
                            letters = letters,
                            OR = OR,
                            as.prob = as.prob,
                            ...)
)
#' @rdname letterFrequency
#' @export
setMethod("letterFrequency", "MaskedModString",
          function(x,
                   letters,
                   OR = "|",
                   as.prob = FALSE)
            letterFrequency(as(x,
                               "ModStringViews"),
                            letters = letters,
                            OR = OR,
                            as.prob = as.prob,
                            collapse = TRUE)
)

# derived from Biostrings/R/letterFrequency.R ----------------------------------

# oligonucleotideFrequency works only for the base codes
# nucleotideFrequencyAt works only for the base codes


# derived from Biostrings/R/letterFrequency.R ----------------------------------
#' @rdname letterFrequency
#' @export
setMethod("consensusMatrix", "ModStringSet",
          function(x,
                   as.prob = FALSE,
                   shift = 0L,
                   width = NULL,
                   baseOnly = FALSE){
            if (!assertive::is_a_bool(as.prob)){
              stop("'as.prob' must be TRUE or FALSE")
            }
            if (!is.integer(shift)){
              shift <- as.integer(shift)
            }
            if (length(x) != 0 && length(shift) > length(x)){
              stop("'shift' has more elements than 'x'")
            }
            if (!is.null(width)) {
              if (!assertive::is_a_number(width) || width < 0L){
                stop("'width' must be NULL or a single non-negative integer")
              }
              if (!is.integer(width)){
                width <- as.integer(width)
              }
            }
            codes <- modscodes(seqtype(x),
                               baseOnly = baseOnly)
            codes[] <- as.integer(unlist(lapply(names(codes),charToRaw)))
            ans <- .Call2("XStringSet_consensus_matrix",
                          x,
                          shift,
                          width,
                          baseOnly,
                          codes,
                          PACKAGE="Biostrings")
            rownames(ans) <- names(modscodes(seqtype(x),
                                             lettersOnly = TRUE))
            ans <- ans[rowSums(ans) > 0, , drop=FALSE]
            if (as.prob) {
              col_sums <- colSums(ans)
              col_sums[col_sums == 0] <- 1  # to avoid division by 0
              ans <- ans / rep(col_sums, each = nrow(ans))
            }
            ans
          }
)

#' @rdname letterFrequency
#' @export
setMethod("consensusString", "ModDNAStringSet",
          function(x,
                   threshold = 0.25,
                   shift = 0L,
                   width = NULL)
            consensusString(consensusMatrix(x,
                                            as.prob = TRUE,
                                            shift = shift,
                                            width = width),
                            ambiguityMap = "?",
                            threshold = threshold)
)

#' @rdname letterFrequency
#' @export
setMethod("consensusString", "ModRNAStringSet",
          function(x,
                   threshold = 0.25,
                   shift = 0L,
                   width = NULL)
            consensusString(consensusMatrix(x,
                                            as.prob = TRUE,
                                            shift = shift,
                                            width = width),
                            ambiguityMap = "?",
                            threshold = threshold)
)

#' @rdname letterFrequency
#' @export
setMethod("consensusString", "ModStringViews",
          function(x,
                   threshold,
                   shift = 0L,
                   width = NULL){
            x <- as(x, "ModStringSet")
            callGeneric()
          }
)
