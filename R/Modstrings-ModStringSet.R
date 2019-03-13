#' @include Modstrings.R
#' @include Modstrings-ModString.R
NULL

# These functions need to be here to access the modified functions of
# - .oneSeqToModStringSet
# - .charToModStringSet
# - format (this dispatches to internal, which screws up the encoding)

#' @name ModStringSet
#' @aliases ModDNAStringSet ModRNAStringSet ModStringSet,ANY-method
#' ModStringSet,AsIs-method ModStringSet,character-method
#' ModStringSet,factor-method ModStringSet,ModString-method
#' ModStringSet,ModStringSet-method ModStringSet,list-method
#' ModStringSet,missing-method as.character,ModStringSet-method
#' show,ModStringSet-method ==,ModStringSet,ModStringSet-method 
#' ==,ModStringSet,XStringSet-method ==,XStringSet,ModStringSet-method
#' 
#' @title ModStringSet objects
#'  
#' @description 
#' The \code{ModStringSet} class is a container for storing a set of 
#' \code{\link{ModString}} objects. It follows the same principles as the 
#' other \code{\link[Biostrings:XStringSet-class]{XStringSet}} objects.
#' 
#' As usual the \code{ModStringSet} containers derive directly from the 
#' \code{\link[Biostrings:XStringSet-class]{XStringSet}} virtual class.
#' 
#' The \code{ModStringSet} class is in itself a virtual class with two types of
#' derivates:
#' \itemize{
#'  \item \code{ModDNAStringSet} 
#'  \item \code{ModRNAStringSet} 
#' }
#' Each class can only be converted to its parent \code{DNAStringSet} or 
#' \code{RNAStringSet}. The modified nucleotides will be converted to their
#' original nucleotides.
#' 
#' Please note, that due to encoding issues not all modifications can be
#' instanciated directly from the console. The vignette contains
#' a comphrensive explanation and examples for working around the problem.
#' 
#' @param x Either a character vector (with no NAs), or an ModString, 
#' ModStringSet or ModStringViews object.
#' @param start,end,width Either NA, a single integer, or an integer vector of 
#' the same length as x specifying how x should be "narrowed" (see ?narrow for 
#' the details).
#' @param use.names TRUE or FALSE. Should names be preserved?
#' 
#' @return a \code{ModStringSet} object.
#' 
#' @examples
#' # Constructing ModDNAStringSet containing an m6A
#' m1 <- ModDNAStringSet(c("AGCT`","AGCT`"))
#' m1
#' 
#' # converting to DNAStringSet
# d1 <- DNAStringSet(m1)
# d2 <- as(m1,"DNAStringSet")
# stopifnot(d1 == d2)
#' 
#' # Constructing ModRNAStringSet containing an m6A
#' m2 <- ModRNAStringSet(c("AGCU`","AGCU`"))
#' m2
NULL


setClass("ModStringSet",contains = c("VIRTUAL","XStringSet"))
#' @rdname ModStringSet
#' @export
setClass("ModDNAStringSet",
         contains = "ModStringSet",
         representation(),
         prototype(
           elementType = "ModDNAString"
         )
)

#' @rdname ModStringSet
#' @export
setClass("ModRNAStringSet",
         contains = "ModStringSet",
         representation(),
         prototype(
           elementType = "ModRNAString"
         )
)

# derived from Biostrings/R/XStringSet-class.R ---------------------------------

#' @export
setReplaceMethod(
  "seqtype", "ModStringSet",
  function(x, value)
  {
    ans_class <- paste0(value, "StringSet")
    if(is(x,ans_class)){
      return(x)
    }
    ans_seq <- .call_new_CHARACTER_from_XStringSet(x)
    ans_seq <- unlist(
      lapply(ans_seq,
             function(a){
               .convert_one_byte_codes_to_originating_base(
                 a,
                 modscodec(seqtype(x)))
             }))
    if (!is.null(names(x))) {
      names(ans_seq) <- names(x)
    }
    do.call(ans_class,list(ans_seq))
  }
)


# derived from Biostrings/R/XStringSet-class.R ---------------------------------

setMethod("make_XStringSet_from_strings", "ModStringSet",
          function(x0, strings, start, end, width)
          {
            codec <- modscodec(seqtype(x0))
            strings <- vapply(strings,
                              function(string)
                                .convert_letters_to_one_byte_codes(string, codec),
                              character(1),
                              USE.NAMES=FALSE)
            callNextMethod()
          }
)

# derived from Biostrings/R/XStringSet-class.R ---------------------------------
# Constructor

#' @rdname ModStringSet
#' @export
ModDNAStringSet <- function(x = character(), start = NA, end = NA, width = NA,
                            use.names = TRUE){
  XStringSet("ModDNA", x, start = start, end = end, width = width,
               use.names = use.names)
  
}
#' @rdname ModStringSet
#' @export
ModRNAStringSet <- function(x = character(), start = NA, end = NA, width = NA,
                            use.names=TRUE){
  XStringSet("ModRNA", x, start = start, end = end, width = width,
              use.names = use.names)
}

# derived from Biostrings/R/XStringSet-class.R ---------------------------------
# Coercion

#' @rdname Modstrings-internals
#' @export
setMethod(
  "XStringSet",
  signature = "ModStringSet",
  function(seqtype, x, start = NA, end = NA, width = NA, use.names = TRUE)
  {
    ans <- narrow(x, start = start, end = end, width = width, 
                  use.names = use.names)
    ans_class <- paste0(seqtype, "StringSet")
    if(is(ans,ans_class)){
      return(ans)
    }
    # convert over "base" classes to convert T/U
    seqtype(ans) <- gsub("Mod","",seqtype(ans))
    seqtype(ans) <- gsub("Mod","",seqtype)
    seqtype(ans) <- seqtype
    ans
  }
)

#' @export
setAs("ANY", "ModDNAStringSet", function(from) ModDNAStringSet(from))
#' @export
setAs("ANY", "ModRNAStringSet", function(from) ModRNAStringSet(from))
#' @export
setMethod(
  "as.character", "ModStringSet",
  function(x, use.names=TRUE)
  {
    ans <- callNextMethod()
    ans <- unlist(lapply(ans,
                         function(a){
                           .convert_one_byte_codes_to_letters(
                             a,
                             modscodec(seqtype(x)))
                         }))
    ans
  }
)

# derived from Biostrings/R/XStringSet-class.R ---------------------------------
# show

.namesW <- 20

.format_utf8 <- function(x,
                         width){
  missingNChar <- width - nchar(x)
  paste(c(x, rep(" ",missingNChar)), collapse = "")
}

.ModStringSet.show_frame_line <- function(x,
                                          i,
                                          iW,
                                          widthW){
  width <- nchar(x)[i]
  snippetWidth <- getOption("width") - 2 - iW - widthW
  if (!is.null(names(x))){
    snippetWidth <- snippetWidth - .namesW - 1
  }
  seq_snippet <- .toSeqSnippet(x[[i]], snippetWidth)
  if (!is.null(names(x))){
    seq_snippet <- .format_utf8(seq_snippet, width = snippetWidth)
  }
  cat(format(paste0("[", i,"]"), width = iW, justify = "right"), " ",
      format(width, width = widthW, justify = "right"), " ",
      seq_snippet,
      sep = "")
  if (!is.null(names(x))) {
    snippet_name <- names(x)[i]
    if (is.na(snippet_name)){
      snippet_name <- "<NA>"
    } else if (nchar(snippet_name) > .namesW) {
      snippet_name <- paste0(substr(snippet_name, 1, .namesW-3), "...")
    }
    cat(" ", snippet_name, sep = "")
  }
  cat("\n")
}

### 'half_nrow' must be >= 1
.ModStringSet.show_frame <- function(x,
                                     half_nrow = 5L){
  if (is.null(head_nrow <- getOption("showHeadLines"))){
    head_nrow <- half_nrow
  }
  if (is.null(tail_nrow <- getOption("showTailLines"))){
    tail_nrow <- half_nrow
  }
  lx <- length(x)
  iW <- nchar(as.character(lx)) + 2 # 2 for the brackets
  ncharMax <- max(nchar(x))
  widthW <- max(nchar(ncharMax), nchar("width"))
  .XStringSet.show_frame_header(iW, widthW, !is.null(names(x)))
  if (lx < (2*half_nrow+1L) | (lx < (head_nrow+tail_nrow+1L))) {
    for (i in seq_len(lx)){
      .ModStringSet.show_frame_line(x, i, iW, widthW)
    }
  } else {
    if (head_nrow > 0){
      for (i in seq_len(head_nrow)){
        .ModStringSet.show_frame_line(x, i, iW, widthW)
      }
    }
    cat(format("...", width = iW, justify = "right"),
        format("...", width = widthW, justify = "right"),
        "...\n")
    if (tail_nrow > 0){
      for (i in (lx-tail_nrow+1L):lx){
        .ModStringSet.show_frame_line(x, i, iW, widthW)
      }
    }
  }
}

#' @export
setMethod(
  "show", "ModStringSet",
  function(object)
  {
    cat("  A ", class(object), " instance of length ",
        length(object), "\n", sep = "")
    if (length(object) != 0){
      .ModStringSet.show_frame(object)
    }
  }
)

# Comparison -------------------------------------------------------------------

.compare_ModStringSet <- function(e1,
                               e2){
  if (!comparable_seqtypes(seqtype(e1), seqtype(e2))) {
    class1 <- class(e1)
    class2 <- class(e2)
    stop("comparison between a \"", class1, "\" instance ",
         "and a \"", class2, "\" instance ",
         "is not supported")
  }
  if(!is(e1,"ModStringSet")){
    e1 <- BStringSet(e1)
  }
  if(!is(e2,"ModStringSet")){
    e2 <- BStringSet(e2)
  }
  pcompare(e1, e2) == 0L
}

#' @export
setMethod("==", signature(e1 = "ModStringSet", e2 = "ModStringSet"),
          function(e1, e2) .compare_ModStringSet(e1, e2)
)
#' @export
setMethod("==", signature(e1 = "ModStringSet", e2 = "XStringSet"),
          function(e1, e2) .compare_ModStringSet(e1, e2)
)
#' @export
setMethod("==", signature(e1 = "XStringSet", e2 = "ModStringSet"),
          function(e1, e2) .compare_ModStringSet(e1, e2)
)
