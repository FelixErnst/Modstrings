#' @include Modstrings.R
NULL

# return the nc independent whether shortName or nomenclature was used
.get_nc_ident <- function(mod,seqtype,nc_type){
  if(nc_type == "short"){
    modNames <- names(modsshortnames(seqtype))
    f <- match(mod,modNames)
    return(names(modsnomenclature(seqtype))[f])
  }
  if(nc_type == "nc"){
    return(mod)
  }
  stop("Something went wrong.",
       call. = FALSE)
}

# normalize the modification type using the nomenclature type
.norm_seqtype_modtype <- function(mod,
                                  seqtype,
                                  nc_type){
  modNames <- switch(nc_type,
                     "short" = modsshortnames(seqtype),
                     "nc" = modsnomenclature(seqtype))
  modValues <- modNames[match(mod,names(modNames))]
  modValues <- modValues[!is.na(modValues)]
  if(length(modValues) != length(mod)){
    stop("No modification for the identifiers '",
         paste(mod[!(mod %in% names(modNames))],
               collapse = "','"),
         "' for a '",paste0(seqtype, "String"),"' found.",
         call. = FALSE)
  }
  modValues
}

#' @rdname modifyNucleotides
#' 
#' @title Modifying nucleotides in a nucleotide sequence (or set of sequences) 
#' at specified locations
#' 
#' @description 
#' \code{\link{modifyNucleotides}} modifies a nucleotide in a sequence (or set 
#' of sequences) based on the type of modification provided. It checks for the 
#' identity of the base nucleotide to be 
#' 
#' @param x a \code{\link{ModString}} or \code{\link{ModStringSet}} object
#' @param at the location where the modification should be made.
#' 
#' The same input as in the original \code{\link[Biostrings]{replaceLetterAt}}
#' are expected:
#'
#' If x is a \code{\link{ModString}} object, then at is typically an integer
#' vector with no NAs but a logical vector or Rle object is valid too. Locations
#' can be repeated and in this case the last replacement to occur at a given
#' location prevails.
#'
#' If x is a rectangular \code{\link{ModStringSet}} object, then \code{at} must
#' be a matrix of logicals with the same dimensions as x. If the
#' \code{\link{ModStringSet}} is not rectangular, \code{at} must be a list of
#' logical vectors.
#' 
#' @param mod The modification short name or nomenclature
#' 
#' If x is a \code{\link{ModString}} object, then letter must be a 
#' \code{\link{ModString}} object or a character vector (with no NAs) with a 
#' total number of letters (sum(nchar(letter))) equal to the number of locations
#' specified in at.
#' 
#' If x is a rectangular \code{\link{ModStringSet}} object, then letter must be
#' a \code{\link{ModStringSet}} object or a character vector of the same length
#' as x. In addition, the number of letters in each element of letter must match
#' the number of locations specified in the corresponding row of at
#' (all(width(letter) == rowSums(at)))
#' 
#' @param nc.type the type of nomenclature to be used. Either "short" or "nc".
#' "Short" for m3C would be "m3C", "nc" for m3C would be "3C". (
#' \code{default = "short"})
#' @param verbose See \code{\link[Biostrings]{replaceLetterAt}}.
#' 
#' @return the input \code{\link{ModString}} or \code{\link{ModStringSet}}
#' object with the changes applied
#' 
#' @export
#' 
#' @examples
#' # modify nucleotides in a ModDNAString 
#' seq <- ModDNAString("AGTC")
#' mseq1 <- modifyNucleotides(seq,c(1,2,4),c("1mA","7mG","3mC"))
#' # This fails since m7G requires a G at the selected position in the sequence
#' \dontrun{
#' mseq <- modifyNucleotides(seq,c(3),c("7mG"))
#' }
#' 
#' # modify nucleotides in a ModRNAString 
#' seq <- ModRNAString("AGUC")
#' mseq1 <- modifyNucleotides(seq,c(1,2,4),c("m1A","m7G","m3C"))
#' # This fails since m7G requires a G at the selected position in the sequence
#' \dontrun{
#' mseq <- modifyNucleotides(seq,c(3),c("m7G"))
#' }
setMethod(
  "modifyNucleotides",
  signature = "ModString",
  definition = function(x, at, mod, nc.type =  c("short","nc"), verbose = FALSE)
  {
    .check_verbose(verbose)
    nc.type <- match.arg(nc.type)
    at <- .check_replace_pos_ModString(x,at)
    assertive::assert_all_are_non_empty_character(mod)
    if(length(at) != length(mod)){
      stop("lengths of 'at' and 'mod' need to be equal.",
           call. = FALSE)
    }
    modValues <- .norm_seqtype_modtype(mod,seqtype(x),nc.type)
    codec <- modscodec(seqtype(x))
    f <- values(codec)[match(modValues, values(codec))]
    current_letter <- unlist(lapply(at,
                                    function(i){
                                      as.character(as(
                                        subseq(x,i,i),
                                        gsub("Mod","",class(x))))
                                    }))
    if(any(originatingBase(codec)[f] != current_letter)){
      mismatch <- which(originatingBase(codec)[f] != current_letter)
      stop("Modification type does not match the originating base:",
           paste("\n",
                 current_letter[mismatch],
                 "!=",
                 originatingBase(codec)[f[mismatch]],
                 " for ",
                 mod[mismatch]),
           call. = FALSE)
    }
    letter <- letters(codec)[f]
    letter <- vapply(letter,
                     .convert_letters_to_one_byte_codes,
                     character(1),
                     modscodec(seqtype(x)))
    for(i in seq_along(at)){
      x <- .call_XString_replace_letter_at(x, at[i], letter[[i]], verbose)
    }
    x
  }
)

#' @rdname modifyNucleotides
#' @export
setMethod(
  "modifyNucleotides",
  signature = "ModStringSet",
  definition = function(x, at, mod, nc.type = c("short","nc"), verbose = FALSE)
  {
    .check_verbose(verbose)
    nc.type <- match.arg(nc.type)
    if (length(x) == 0L){
      stop("'x' has no element")
    }
    .norm_replace_pos_ModStringSet(x,at)
    if (is(mod, "ModStringSet")) {
      mod <- lapply(mod,
                    function(m){
                      unname(as.character(split(m,seq_len(length(m)))))
                    })
    }
    if (!is.list(mod)){
      stop("'mod' must be a ModStringSet object or a list of character vectors")
    }
    .check_letter_ModStringSet(x,at,mod)
    unlisted_x <- unlist(x, use.names = FALSE)
    if(is.list(at)){
      at <- unlist(at)
    } else {
      at <- as.vector(t(at))
    }
    unlisted_ans <- modifyNucleotides(unlisted_x,
                                      at,
                                      unlist(mod),
                                      nc.type = nc.type,
                                      verbose = verbose)
    relist(unlisted_ans, x)
  }
)

#' @rdname modifyNucleotides
#' @export
setMethod(
  "modifyNucleotides",
  signature = "DNAString",
  definition = function(x, at, mod, nc.type = c("short","nc"), verbose = FALSE)
  {
    modifyNucleotides(as(x,"ModDNAString"), at, mod, nc.type = nc.type,
                      verbose = verbose)
  }
)
#' @rdname modifyNucleotides
#' @export
setMethod(
  "modifyNucleotides",
  signature = "RNAString",
  definition = function(x, at, mod, nc.type = c("short","nc"), verbose = FALSE)
  {
    modifyNucleotides(as(x,"ModRNAString"), at, mod, nc.type = nc.type,
                      verbose = verbose)
  }
)
#' @rdname modifyNucleotides
#' @export
setMethod(
  "modifyNucleotides",
  signature = "DNAStringSet",
  definition = function(x, at, mod, nc.type = c("short","nc"), verbose = FALSE)
  {
    modifyNucleotides(as(x,"ModDNAStringSet"), at, mod, nc.type = nc.type,
                      verbose = verbose)
  }
)
#' @rdname modifyNucleotides
#' @export
setMethod(
  "modifyNucleotides",
  signature = "RNAStringSet",
  definition = function(x, at, mod, nc.type = c("short","nc"), verbose = FALSE)
  {
    modifyNucleotides(as(x,"ModRNAStringSet"), at, mod, nc.type = nc.type,
                      verbose = verbose)
  }
)
