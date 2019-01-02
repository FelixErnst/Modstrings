#' @include Modstrings.R
#' @include Modstrings-ModStringSet.R
NULL

#' @name ModStringSetList
#' @aliases ModDNAStringSetList ModRNAStringSetList
#' 
#' @title ModStringSetList
#' 
#' @description 
#' title
#' 
#' @param ... \code{\link{ModStringSet}} objects of one type. 
#' @param use.names \code{TRUE}(default) or \code{FALSE}: Whether names of the 
#' input \code{ModStringSet} objects should be stored and used as the element
#' names in the \code{ModStringSetList}.
#' 
#' @return a \code{ModStringSetList} object.
#' 
#' @examples 
#' mrseq <- c("ACGU7","ACGU7","ACGU7","ACGU7")
#' 
#' # Example: contruction of ModStringSetlist from ModString objects
#' mr <- ModRNAString("ACGU7")
#' mrs <- ModRNAStringSet(list(mr,mr,mr,mr))
#' mrsl <- ModRNAStringSetList(mrs,mrs)
#' 
#' # Example: construction of ModStringSetlist from mixed sources
#' mrsl2 <- ModRNAStringSetList(mrs,mrseq)
#' 
NULL

# These functions need to be here to access the modified functions of
# - .new_ModStringSetList_from_list
# - .new_ModStringSetList_from_List

# derived from Biostrings/R/XStringSetList-class.R -----------------------------

setClass("ModStringSetList", contains = "XStringSetList")
#' @rdname ModStringSetList
#' @export
setClass("ModDNAStringSetList",
         contains = "ModStringSetList",
         representation(
           unlistData = "ModDNAStringSet"
         ),
         prototype(
           elementType = "ModDNAStringSet"
         )
)
#' @rdname ModStringSetList
#' @export
setClass("ModRNAStringSetList",
         contains = "ModStringSetList",
         representation(
           unlistData = "ModRNAStringSet"
         ),
         prototype(
           elementType = "ModRNAStringSet"
         )
)

# derived from Biostrings/R/XStringSetList-class.R -----------------------------

.new_ModStringSetList_from_list <- function(seqtype,
                                            x){
  x_eltNROWS <- S4Vectors::elementNROWS(x)
  empty_idx <- which(x_eltNROWS == 0L)
  if (length(empty_idx) != 0L) {
    y <- x[-empty_idx]
  } else {
    y <- x
  }
  unlisted_y <- unlist(y, use.names=FALSE, recursive=FALSE)
  if (!is.list(unlisted_y) && length(unlisted_y) == sum(x_eltNROWS)) {
    unlisted_ans <- ModStringSet(seqtype,
                                 unlisted_y)
  } else {
    ## In that case 'length(unlisted_y)' should be < 'sum(x_eltNROWS)'
    ## which means unlist() was not able to fully unlist 'y'. So let's
    ## try to turn each list element into an XStringSet object and then
    ## combine them together. This is of course much slower than if
    ## unlist() had succeeded.
    y <- lapply(unname(y),
                ModStringSet,
                seqtype = seqtype)
    unlisted_ans <- do.call(c, y)
  }
  relist(unlisted_ans, x)
}

.new_ModStringSetList_from_List <- function(seqtype,
                                            x){
  unlisted_x <- unlist(x,
                       use.names = FALSE)
  unlisted_ans <- ModStringSet(seqtype,
                               unlisted_x)
  ans <- relist(unlisted_ans, x)
  ## relist() puts the names back but not the metadata columns.
  mcols(ans) <- mcols(x)
  ans
}

ModStringSetList <- function(seqtype,
                             ...,
                             use.names = TRUE){
  if (!assertive::is_a_bool(use.names)){
    stop("'use.names' must be TRUE or FALSE")
  }
  x <- list(...)
  if (length(x) == 1L) {
    x1 <- x[[1L]]
    if (is.list(x1) || (is(x1, "List") && !is(x1, "ModStringSet"))) {
      x <- x1
      if (is(x, "List")) {
        if (!use.names)
          names(x) <- NULL
        return(.new_ModStringSetList_from_List(seqtype,
                                               x))
      }
    }
  }
  if (!use.names)
    names(x) <- NULL
  .new_ModStringSetList_from_list(seqtype,
                                  x)
}

# derived from Biostrings/R/XStringSetList-class.R -----------------------------
# Constructor

#' @rdname ModStringSetList
#' @export
ModDNAStringSetList <- function(...,
                                use.names = TRUE){
  ModStringSetList("ModDNA",
                   ...,
                   use.names = use.names)
}
#' @rdname ModStringSetList
#' @export
ModRNAStringSetList <- function(...,
                                use.names = TRUE){
  ModStringSetList("ModRNA",
                   ...,
                   use.names = use.names)
}
