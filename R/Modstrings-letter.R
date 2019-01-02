#' @include Modstrings.R
#' @include Modstrings-ModString.R
#' @include Modstrings-ModStringViews.R
NULL

#' @name letter
#' 
#' @title Subsetting a ModString to a single letter
#' 
#' @description 
#' As the \code{\link[Biostrings:letter]{Biostrings}} functions, extracts a
#' substring from a \code{\link{ModString}} object by returning individual
#' letters by their position as a character vector.
#'
#' @param x ModString 
#' @param i the position of the letter to be returned 
#'
#' @return a character vector
#' @export
#'
#' @examples
#' seq <- ModDNAString("AGCT/")
#' letter(seq, 5)
NULL

# These functions need to be here to access the modified functions of
# - ModString.read

# derived from Biostrings/R/letter.R -------------------------------------------
### Return a character vector of length 1.
#' @rdname letter
#' @export
setMethod("letter", "ModString",
          function(x, i)
          {
            if (!is.numeric(i) || any(is.na(i)))
              stop("'i' must be an NA-free numeric vector")
            if (!all(i >= 1) || !all(i <= x@length))
              stop("subscript out of bounds")
            ModString.read(x, i)
          }
)
### Return a character vector of the same length as 'x'.
#' @rdname letter
#' @export
setMethod("letter", "ModStringViews",
          function(x, i){
            if (!is.numeric(i) || any(is.na(i))){
              stop("'i' must be an NA-free numeric vector")
            }
            if (length(x) == 0){
              return(character(0))
            }
            imax <- min(nchar(x))
            if (!all(i >= 1) || !all(i <= imax)){
              stop("subscript out of bounds")
            }
            vapply(seq_len(length(x)), 
                   function(n) {
                     ModString.read(x[[n]], i)
                   },
                   character(1))
          }
)
