#' @include Modstrings.R
NULL

# normalize the modification type using the nomenclature type
.norm_seqtype_modtype <- function(mod,
                                  seqtype,
                                  nc_type,
                                  class){
  modNames <- switch(nc_type,
                     "short" = modsshortnames(seqtype),
                     "nc" = modsnomenclature(seqtype))
  modValues <- modNames[match(mod,names(modNames))]
  modValues <- modValues[!is.na(modValues)]
  if(length(modValues) != length(mod)){
    stop("No modification for the identifiers '",
         paste(mod[!(mod %in% names(modNames))],
               collapse = "','"),
         "' for a '",class,"' found.",
         call. = FALSE)
  }
  modValues
}

.check_stop.on.error <- function(stop.on.error)
{
  if (!.is_a_bool(stop.on.error)){
    stop("'stop.on.error' must be TRUE or FALSE",
         call. = FALSE)
  }
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
#' If \code{x} is a \code{\link{ModString}} object, then letter must be a 
#' \code{\link{ModString}} object or a character vector (with no \code{NA}) with
#' a total number of letters \code{(sum(nchar(letter)))} equal to the number of
#' locations specified in at.
#' 
#' If \code{x} is a rectangular \code{\link{ModStringSet}} object, then letter
#' must be a \code{\link{ModStringSet}} object, a list of character vectors or a
#' \code{\link[IRanges:AtomicList-class]{CharacterList}} of the same length as
#' x. In addition, the number of letters in each element of letter must match
#' the number of locations specified in the corresponding row of at
#' \code{(all(width(letter) == rowSums(at)))}.
#' 
#' @param nc.type the type of nomenclature to be used. Either "short" or "nc".
#' "Short" for m3C would be "m3C", "nc" for m3C would be "3C". (
#' \code{default = "short"})
#' @param stop.on.error For \code{combineIntoModstrings}: \code{TRUE}(default)
#' or \code{FALSE}: Should an error be raised upon encounter of incompatible 
#' positions?
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
#' seq
#' 
#' mseq1 <- modifyNucleotides(seq,c(1,2,4),c("1mA","7mG","3mC"))
#' mseq1
#' 
#' # This fails since m7G requires a G at the selected position in the sequence
#' \dontrun{
#' mseq <- modifyNucleotides(seq,c(3),c("7mG"))
#' }
#' 
#' # modify nucleotides in a ModRNAString 
#' seq <- ModRNAString("AGUC")
#' seq
#' 
#' mseq1 <- modifyNucleotides(seq,c(1,2,4),c("m1A","m7G","m3C"))
#' mseq1
#' 
#' # This fails since m7G requires a G at the selected position in the sequence
#' \dontrun{
#' mseq <- modifyNucleotides(seq,c(3),c("m7G"))
#' }
setMethod(
  "modifyNucleotides",
  signature = "ModString",
  definition = function(x, at, mod, nc.type =  c("short","nc"),
                        stop.on.error = TRUE, verbose = FALSE)
  {
    .check_verbose(verbose)
    .check_stop.on.error(stop.on.error)
    nc.type <- match.arg(nc.type)
    at <- .norm_replace_pos_ModString(x,at)
    if(!.is_non_empty_character(as.character(unlist(mod)))){
      stop("'mod' must all be non empty characters.", call. = FALSE)
    }
    if(length(at) != length(mod)){
      stop("lengths of 'at' and 'mod' need to be equal.",
           call. = FALSE)
    }
    # check if originating base matches the modification
    modValues <- .norm_seqtype_modtype(mod, seqtype(x), nc.type, class(x))
    codec <- modscodec(seqtype(x))
    f <- values(codec)[match(modValues, values(codec))]
    current_letter <- lapply(at,
                             function(i){
                               subseq(x,i,i)
                             })
    current_letter <- vapply(current_letter,as.character,character(1))
    class <- paste0(class(x),"Set")
    current_letter <- as(do.call(class, list(current_letter)),
                         gsub("Mod","",class))
    current_letter <- as.character(current_letter)
    if(any(originatingBase(codec)[f] != current_letter)){
      mismatch <- originatingBase(codec)[f] != current_letter
      n <- min(5L, sum(mismatch))
      msg <- paste0("Modification type does not match the originating base:",
                    paste("\n",
                          current_letter[mismatch][seq_len(n)],
                          "!=",
                          originatingBase(codec)[f[mismatch][seq_len(n)]],
                          " for ",
                          mod[mismatch][seq_len(n)],
                          collapse = ""))
      if(sum(mismatch) > 5L){
        msg <- paste0(msg,"\nand more ...")
      }
      if(stop.on.error){
        stop(msg, call. = FALSE)
      } else {
        warning(paste0(msg,"\nSkipping these position(s) ..."), call. = FALSE)
      }
    } else {
      mismatch <- rep(FALSE,length(f))
    }
    #
    letter <- letters(codec)[f]
    letter <- vapply(letter,
                     .convert_letters_to_one_byte_codes,
                     character(1),
                     modscodec(seqtype(x)))
    for(i in seq_along(at[!mismatch])){
      x <- .call_XString_replace_letter_at(x, at[!mismatch][i],
                                           letter[!mismatch][[i]], verbose)
    }
    x
  }
)

#' @rdname modifyNucleotides
#' @export
setMethod(
  "modifyNucleotides",
  signature = "ModStringSet",
  definition = function(x, at, mod, nc.type = c("short","nc"), 
                        stop.on.error = TRUE, verbose = FALSE)
  {
    .check_verbose(verbose)
    .check_stop.on.error(stop.on.error)
    nc.type <- match.arg(nc.type)
    if (length(x) == 0L){
      stop("'x' has no element")
    }
    .check_replace_pos_ModStringSet(x,at)
    if (is(mod, "ModStringSet")) {
      tmp <- separate(mod, nc.type = nc.type)
      mod <- split(unname(mcols(tmp)$mod), seqnames(tmp))
    }
    if (!is.list(mod) && !is(mod,"CharacterList")){
      stop("'mod' must be a ModStringSet object, a list of character vectors ",
           "or a CharacterList of the same length as 'x'",
           call. = FALSE)
    }
    # a preamptive check - return value is not used
    .norm_seqtype_modtype(unlist(mod), seqtype(x), nc.type, class(x))
    #
    .check_letter_ModStringSet(x,at,mod)
    unlisted_x <- unlist(x, use.names = FALSE)
    if(is.list(at)){
      at <- unlist(IRanges::LogicalList(at)) # apparently a bit more efficient
    } else {
      at <- as.vector(t(at))
    }
    unlisted_ans <- modifyNucleotides(unlisted_x,
                                      at,
                                      unlist(mod),
                                      nc.type = nc.type,
                                      stop.on.error = stop.on.error, 
                                      verbose = verbose)
    relist(unlisted_ans, x)
  }
)

#' @rdname modifyNucleotides
#' @export
setMethod(
  "modifyNucleotides",
  signature = "DNAString",
  definition = function(x, at, mod, nc.type = c("short","nc"), 
                        stop.on.error = TRUE, verbose = FALSE)
  {
    modifyNucleotides(as(x,"ModDNAString"), at, mod, nc.type = nc.type,
                      stop.on.error = stop.on.error, verbose = verbose)
  }
)
#' @rdname modifyNucleotides
#' @export
setMethod(
  "modifyNucleotides",
  signature = "RNAString",
  definition = function(x, at, mod, nc.type = c("short","nc"),
                        stop.on.error = TRUE, verbose = FALSE)
  {
    modifyNucleotides(as(x,"ModRNAString"), at, mod, nc.type = nc.type,
                      stop.on.error = stop.on.error, verbose = verbose)
  }
)
#' @rdname modifyNucleotides
#' @export
setMethod(
  "modifyNucleotides",
  signature = "DNAStringSet",
  definition = function(x, at, mod, nc.type = c("short","nc"),
                        stop.on.error = TRUE, verbose = FALSE)
  {
    modifyNucleotides(as(x,"ModDNAStringSet"), at, mod, nc.type = nc.type,
                      stop.on.error = stop.on.error, verbose = verbose)
  }
)
#' @rdname modifyNucleotides
#' @export
setMethod(
  "modifyNucleotides",
  signature = "RNAStringSet",
  definition = function(x, at, mod, nc.type = c("short","nc"),
                        stop.on.error = TRUE, verbose = FALSE)
  {
    modifyNucleotides(as(x,"ModRNAStringSet"), at, mod, nc.type = nc.type,
                      stop.on.error = stop.on.error, verbose = verbose)
  }
)
