#' @include Modstrings.R
NULL

#' @name replaceAt
#' @aliases extractAt extractAt,ModString-method extractAt,ModStringSet-method
#' replaceAt,ModString-method replaceAt,ModStringSet-method
#' 
#' @title Extract/replace a substring from/in a ModString or ModStringSet
#' 
#' @description 
#' \code{extractAt} extracts multiple subsequences from a 
#' \code{\link{ModString}} object x, or from the individual sequences of 
#' \code{\link{ModStringSet}} object x, at the ranges of positions specified 
#' thru at.
#' 
#' \code{replaceAt} performs multiple subsequence replacements 
#' (a.k.a. substitutions) in \code{\link{ModString}} object x, or in the
#' individual sequences of \code{\link{ModStringSet}} object x,
#' at the ranges of positions specified thru at.
#' 
#' @param x a \code{\link{ModString}} or a \code{\link{ModStringSet}} object.
#' @param at See \code{\link[Biostrings]{replaceAt}}.
#' @param value a \code{character} vector, \code{\link{ModString}} or a 
#' \code{\link{ModStringSet}} object. See \code{\link[Biostrings]{replaceAt}}.
#'
#' @return the modified \code{\link{ModString}} or \code{\link{ModStringSet}} 
#' object
#' 
#' @export
#'
#' @examples
#' seq1 <- ModDNAString("AGCTTTT//TTTCGA")
#' seq2 <- replaceAt(seq1,IRanges::IRanges(7,10),"AGCT")
#' set1 <- ModDNAStringSet(c("AGCTTTT//TTTTCGA","AGCTTT//TTTTTCGA"))
#' set2 <- replaceAt(set1,IRanges::IRanges(7,10),"AGCT")
NULL

# derived from Biostrings/R/replaceAt.R ----------------------------------------

.make_ModStringSet_from_value <- function(value, x_seqtype)
{
  if (!is(value, "ModStringSet")) {
    value_class <- paste0(x_seqtype, "StringSet")
    value <- try(as(value, value_class), silent=TRUE)
    if (is(value, "try-error")){
      stop("failed to coerce 'value' to a ", value_class, " object")
    }
  }
  value
}

.make_ModStringSetList_from_value <- function(value, x_seqtype)
{
  if (is.character(value)) {
    value_class <- paste0(x_seqtype, "StringSet")
    value <- as(value, value_class)
  }
  if (is(value, "ModStringSet")) {
    value <- relist(value, list(seq_along(value)))
  }
  if (is.list(value)){
    value <- IRanges::CharacterList(value)
  }
  if (is(value, "CharacterList")) {
    unlisted_value <- unlist(value, use.names=FALSE)
    unlisted_value <- .make_ModStringSet_from_value(unlisted_value, x_seqtype)
    value <- relist(unlisted_value, value)
  }
  if (seqtype(value) != x_seqtype){
    seqtype(value) <- x_seqtype
  }
  if (!is(value, "ModStringSetList")){
    stop("invalid type of 'value'")
  }
  value
}

.normarg_value1 <- function(value, at, x_seqtype)
{
  value <- .make_ModStringSet_from_value(value, x_seqtype)
  Biostrings:::.V_recycle(value,
                          at,
                          "value",
                          "the number of replacements")
}

### 'at' is assumed to be normalized so it has the length of 'x'.
.normarg_value2 <- function(value, at, x_seqtype)
{
  value <- .make_ModStringSetList_from_value(value, x_seqtype)
  value <- Biostrings:::.V_recycle(value,
                                   at,
                                   "value",
                                   "'length(x)'")
  Biostrings:::.H_recycle(value,
                          at,
                          "value",
                          "at",
                          paste0("after recycling of 'at' and 'value' ",
                                 "to the length of 'x'"))
}


# derived from Biostrings/R/replaceAt.R ----------------------------------------
# extractAt()

#' @rdname replaceAt
#' @export
setMethod("extractAt", "ModString",
          function(x, at)
          {
            at <- Biostrings:::.make_IRanges_from_at(at)
            ## extractList() will check that all the ranges in 'at' are within
            ## the limits of sequence 'x'.
            IRanges::extractList(x, at)
          }
)
#' @rdname replaceAt
#' @export
setMethod("extractAt", "ModStringSet",
          function(x, at)
          {
            at <- Biostrings:::.normarg_at2(at, x)
            at_eltNROWS <- elementNROWS(at)
            x2 <- rep.int(unname(x), at_eltNROWS)
            unlisted_at <- unlist(at, use.names=FALSE)
            unlisted_ans <- subseq(x2,
                                   start = start(unlisted_at),
                                   width = width(unlisted_at))
            names(unlisted_ans) <- names(unlisted_at)
            ans <- relist(unlisted_ans, at)
            names(ans) <- names(x)
            ans
          }
)


# derived from Biostrings/R/replaceAt.R ----------------------------------------
# replaceAt()

#' @rdname replaceAt
#' @export
setMethod("replaceAt", "ModString",
          function(x, at, value = ""){
            if (length(at) == 0L && length(value) == 0L){
              return(x)
            }
            at <- Biostrings:::.normarg_at1(at, x)
            value <- .normarg_value1(value, at, seqtype(x))
            NR <- length(at)  # same as length(value) -- nb of replacements
            if (NR == 0L){
              return(x)
            }
            .Call2("XString_replaceAt",
                   x,
                   at,
                   value,
                   PACKAGE = "Biostrings")
          }
)
#' @rdname replaceAt
#' @export
setMethod("replaceAt", "ModStringSet",
          function(x, at, value=""){
            if (length(at) == 0L && length(value) == 0L){
              return(x)
            }
            at <- Biostrings:::.normarg_at2(at, x)
            value <- .normarg_value2(value, at, seqtype(x))
            ans <- .Call2("XStringSet_replaceAt",
                          x,
                          at,
                          value,
                          PACKAGE = "Biostrings")
            names(ans) <- names(x)
            mcols(ans) <- mcols(x)
            ans
          }
)
