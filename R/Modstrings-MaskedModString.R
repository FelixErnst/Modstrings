#' @include Modstrings.R
#' @include Modstrings-ModStringSet.R
#' @include Modstrings-ModStringViews.R
NULL

#' @name MaskedModString
#' @aliases seqtype,MaskedModString-method
#' 
#' @title MaskedModString objects
#'  
#' @description 
#' The functions are implemented as defined in the Biostrings package. Have
#' a look the \code{\link[Biostrings:MaskedXString-class]{MaskedXString}} class.
#' 
#' @param x a \code{ModString} object.
#' @param start,end not used for \code{MaskedModString} objects.
#' 
#' @return a \code{MaskedModString} object.
#' 
#' @examples
#' # Mask positions
#' mask <- Mask(mask.width=5, start=c(2), width=c(3))
#' mr <- ModRNAString("ACGU7")
#' masks(mr) <- mask
#' 
#' # Invert masks
#' mr <- gaps(mr)
#' 
#' # Drop the mask
#' masks(mr) <- NULL
NULL

setClass("MaskedModString",
         contains =  c("MaskedXString"),
         representation("VIRTUAL",
                        unmasked = "ModString")
)
setClass("MaskedModDNAString",
         contains = "MaskedModString",
         representation(
           unmasked = "ModDNAString"
         )
)
setClass("MaskedModRNAString",
         contains = "MaskedModString",
         representation(
           unmasked = "ModRNAString"
         )
)

# derived from Biostrings/R/MaskedXString.R ------------------------------------
# Accessor-like methods.
# should work from Biostrings

# derived from Biostrings/R/MaskedXString.R ------------------------------------
# Validity


# derived from Biostrings/R/MaskedXString.R ------------------------------------
# "seqtype" and "seqtype<-" methods

#' @rdname MaskedModString
#' @export
setMethod("seqtype", "MaskedModString", function(x) seqtype(unmasked(x)))

# derived from Biostrings/R/MaskedXString.R ------------------------------------
# Coercion 

### From MaskedXString objects to MaskedModString objects.
setAs("MaskedXString", "MaskedModDNAString",
      function(from) {seqtype(from) <- "ModDNA"; from}
)
setAs("MaskedXString", "MaskedModRNAString",
      function(from) {seqtype(from) <- "ModRNA"; from}
)

setAs("ModDNAString", "MaskedModDNAString",
      function(from){
        masks <- new("MaskCollection", width=length(from))
        new("MaskedModDNAString", unmasked=from, masks=masks)
      }
)
setAs("ModRNAString", "MaskedModRNAString",
      function(from){
        masks <- new("MaskCollection", width = length(from))
        new("MaskedModRNAString", unmasked = from, masks = masks)
      }
)

### From MaskedModString objects to ModString objects.
setAs("MaskedModDNAString", "ModDNAString",
      function(from) unmasked(from)
)
setAs("MaskedModRNAString", "ModRNAString",
      function(from) unmasked(from)
)

setMethod("ModString", "MaskedModString",
          function(seqtype, x, start=NA, end=NA, width=NA)
            ModString(seqtype, unmasked(x), start=start, end=end, width=width)
)

### From a MaskedXString object to a MaskCollection object.
setAs("MaskedModString", "MaskCollection",
      function(from) masks(from)
)

### From a MaskedXString object to a NormalIRanges object.
setAs("MaskedModString", "NormalIRanges",
      function(from) as(masks(from), "NormalIRanges")
)

### From a MaskedXString object to an XStringViews object.
setAs("MaskedModString", "ModStringViews",
      function(from)
      {
        views <- gaps(collapse(masks(from)))[[1]]
        unsafe.newModStringViews(unmasked(from), start(views), width(views))
      }
)

setAs("MaskedModString", "Views", function(from) as(from, "ModStringViews"))

### NOT exported.
toModStringViewsOrModString <- function(x){
  x0 <- unmasked(x)
  mask1 <- collapse(masks(x))
  if (S4Vectors::isEmpty(mask1))
    return(x0)
  views <- gaps(mask1)[[1]]
  unsafe.newModStringViews(x0, start(views), width(views))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The transformation methods (endomorphisms) "collapse" and "gaps".
###

#' @rdname MaskedModString
#' @export
setMethod("collapse", "MaskedXString",
          function(x){
            x@masks <- collapse(masks(x))
            x
          }
)

### 'start' and 'end' are ignored.
#' @rdname MaskedModString
#' @export
setMethod("gaps", "MaskedXString",
          function(x, start = NA, end = NA){
            x@masks <- gaps(collapse(masks(x)))
            x
          }
)
