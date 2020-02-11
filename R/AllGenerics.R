#' @include Modstrings.R
NULL

# accessors for ModCodec -------------------------------------------------------

setGeneric(name = "letters",
           signature = "x",
           def = function(x) standardGeneric("letters"))
setGeneric(name = "oneByteCodes",
           signature = "x",
           def = function(x) standardGeneric("oneByteCodes"))
setGeneric(name = "conversion",
           signature = "x",
           def = function(x) standardGeneric("conversion"))
setGeneric(name = "originatingBase",
           signature = "x",
           def = function(x) standardGeneric("originatingBase"))
setGeneric(name = "values",
           signature = "x",
           def = function(x) standardGeneric("values"))
setGeneric(name = "lettersEscaped",
           signature = "x",
           def = function(x) standardGeneric("lettersEscaped"))
setGeneric(name = "oneByteCodesEscaped",
           signature = "x",
           def = function(x) standardGeneric("oneByteCodesEscaped"))
setGeneric(name = "lettersNeedEscape",
           signature = "x",
           def = function(x) standardGeneric("lettersNeedEscape"))
setGeneric(name = "oneByteCodesNeedEscape",
           signature = "x",
           def = function(x) standardGeneric("oneByteCodesNeedEscape"))
setGeneric(name = "additionalInfo",
           signature = "x",
           def = function(x) standardGeneric("additionalInfo"))

# seqtype ----------------------------------------------------------------------

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
           def = function(x, at, mod, nc.type = "short", stop.on.error = TRUE,
                          verbose = FALSE)
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
                          stop.on.error = TRUE, verbose = FALSE, ...)
             standardGeneric("combineIntoModstrings"))

#' @rdname separate
#' @export
setGeneric(name = "removeIncompatibleModifications",
           signature = c("gr", "x"),
           def = function(gr, x, ...)
             standardGeneric("removeIncompatibleModifications"))
