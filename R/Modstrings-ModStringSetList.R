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

# Constructor ------------------------------------------------------------------

#' @rdname ModStringSetList
#' @export
ModDNAStringSetList <- function(..., use.names = TRUE){
  XStringSetList("ModDNA", ..., use.names = use.names)
}
#' @rdname ModStringSetList
#' @export
ModRNAStringSetList <- function(..., use.names = TRUE){
  XStringSetList("ModRNA", ..., use.names = use.names)
}
