#' @include Modstrings.R
#' @include Modstrings-ModStringSet.R
NULL

#' @name ModStringViews
#' @aliases ModStringViews show,ModStringViews-method 
#' ==,ModStringViews,ModStringViews-method ==,ModStringViews,XString-method 
#' ==,XStringViews,ModString-method ModStringSet,ModStringViews-method
#' 
#' @title The ModStringViews class extending the XStringViews class
#' 
#' @description 
#' As the \code{\link[Biostrings:XStringViews-class]{XStringViews}} the 
#' \code{ModStringViews} is the basic container for storing a set of views on 
#' the same sequence (this time a \code{ModString} object).
#' 
#' @details 
#' For the details have a look at the 
#' \code{\link[Biostrings:XStringViews-class]{XStringViews}} class.
#' 
#' @param subject,start,end,width,names
#' See \code{\link[Biostrings:XStringViews-class]{XStringViews}}.
#' 
#' @return a \code{ModStringViews} object.
#' 
#' @examples
#' seq <- ModDNAString("AGC6AGC6")
#' v <- Views(seq, start = 3:1, end = 6:8)
NULL

# derived from Biostrings/R/XStringViews-class.R -------------------------------
setClass("ModStringViews",
         contains="XStringViews",
         representation(
           subject = "ModString"
         )
)

# Constructor

#' @rdname ModStringViews
#' @export
setMethod(
  "Views", "ModString",
  function(subject, start = NULL, end = NULL, width = NULL, 
           names = NULL)
  {
    .new_Views(subject, start = start, end = end, width = width, names = names,
               Class = "ModStringViews")
  }
)

# Coercion

#' @export
setAs("ModStringViews", "ModDNAStringSet", function(from) ModDNAStringSet(from))
#' @export
setAs("ModStringViews", "ModRNAStringSet", function(from) ModRNAStringSet(from))
#' @export
setAs("ModStringSet", "Views", .XStringSetAsViews)
#' @export
setAs("ModStringSet", "ModStringViews", .XStringSetAsViews)
#' @export
setAs("ModStringViews", "ModDNAStringSet", function(from) ModDNAStringSet(from))
#' @export
setAs("ModStringViews", "ModRNAStringSet", function(from) ModRNAStringSet(from))


# Comparison

# These functions need to be here to access the modified functions of
# - comparable_seqtypes

#' @export
setMethod("==", signature(e1 = "ModStringViews", e2 = "ModStringViews"),
          function(e1, e2){
            if (!comparable_seqtypes(seqtype(e1), seqtype(e2))) {
              class1 <- class(subject(e1))
              class2 <- class(subject(e2))
              stop("comparison between XStringViews objects with subjects of ",
                   "class \"", class1, "\" and \"", class2, "\" ",
                   "is not supported")
            }
            .XStringViews.equal(e1, e2)
          }
)

.compare_ModString_Views <- function(e1,
                                     e2){
  if (!comparable_seqtypes(seqtype(e1), seqtype(e2))) {
    class1 <- class(subject(e1))
    class2 <- class(e2)
    stop("comparison between an XStringViews object with a subject of ",
         "class \"", class1, "\" and a \"", class2, "\" instance ",
         "is not supported")
  }
  .XStringViews.equal(e1, as(e2, "Views"))
}
#' @export
setMethod("==", signature(e1 = "ModStringViews", e2 = "XString"), 
          function(e1, e2) .compare_ModString_Views(e1, e2)
)
#' @export
setMethod("==", signature(e1 = "XStringViews", e2 = "ModString"), 
          function(e1, e2) .compare_ModString_Views(e1, e2)
)
