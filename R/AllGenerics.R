#' @include Modstrings.R
NULL

# seqtype -----------------------------------------------------------------------

#' @rdname alphabet
#' @export
setGeneric(name = "shortName",
           signature = "x",
           def = function(x) standardGeneric("shortName"))
#' @rdname alphabet
#' @export
setGeneric(name = "fullName",
           signature = "x",
           def = function(x) standardGeneric("fullName"))
#' @rdname alphabet
#' @export
setGeneric(name = "nomenclature",
           signature = "x",
           def = function(x) standardGeneric("nomenclature"))


# modify -----------------------------------------------------------------------

#' @name modifyNucleotides
#' @export
setGeneric(name = "modifyNucleotides",
           signature = "x",
           def = function(x, at, mod, nc.type = "short", verbose = FALSE)
             standardGeneric("modifyNucleotides"))


# seperate ---------------------------------------------------------------------

#' @name separate
#' @export
setGeneric(name = "separate",
           signature = "x",
           def = function(x, nc.type = "short") standardGeneric("separate"))

#' @rdname separate
#' @export
setGeneric(name = "combineIntoModstrings",
           signature = c("x","gr"),
           def = function(x, gr, with.qualities = FALSE, quality.type = "Phred",
                          verbose = FALSE, ...)
             standardGeneric("combineIntoModstrings"))
