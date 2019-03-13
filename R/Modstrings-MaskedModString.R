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

### From a MaskedXString object to an XStringViews object.
setAs("MaskedModString", "ModStringViews",
      function(from)
      {
        views <- gaps(collapse(masks(from)))[[1]]
        Views(unmasked(from), start = start(views), width = width(views))
      }
)
