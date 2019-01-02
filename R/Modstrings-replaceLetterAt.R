#' @include Modstrings.R
NULL

# derived from Biostrings/R/replaceLetterAt.R ----------------------------------

.check_replace_pos_ModString <- function(x,
                                         at){
  if (is(at, "Rle")){
    at <- as.vector(at)
  }
  if (is.logical(at)) {
    if (length(at) != length(x)){
      stop("when 'at' is a logical sequence, it must have the ",
           "same length as 'x'",
           call. = FALSE)
    }
    at <- which(at)
  } else {
    if (!is.numeric(at)){
      stop("'at' must be a vector of integers",
           call. = FALSE)
    }
    if (!is.integer(at)){
      at <- as.integer(at)
    }
  }
  at
}

.check_replace_pos_ModStringSet <- function(x,at){
  x_width <- width(x)
  if(is.list(at)){
    if(!is.logical(unlist(at))){
      stop("'at' must be a matrix or list of logicals",
           call. = FALSE)
    }
    if (length(at) != length(x) || lengths(at) != x_width){
      stop("'x' and 'at' must have the same dimensions",
           call. = FALSE)
    }
  } else {
    if(!is.logical(at) || !is.matrix(at)){
      stop("'at' must be a matrix or list of logicals",
           call. = FALSE)
    }
    if (nrow(at) != length(x) || ncol(at) != x_width){
      stop("'x' and 'at' must have the same dimensions",
           call. = FALSE)
    }
  }
}

.check_letter_ModStringSet <- function(
  x,
  at,
  letter,
  .xname = assertive::get_name_in_parent(letter))
{
  if (length(letter) != length(x)){
    stop("'x' and '",.xname,"' must have the same length",
         call. = FALSE)
  }
  if(is.list(letter)){
    lengths <- lengths(letter)
  } else {
    lengths <- width(letter)
  }
  if(is.matrix(at)){
    sum <- rowSums(at)
  } else if(is.list(at)){
    sum <- vapply(at,sum,numeric(1))
  } else {
    stop("Something went wrong.")
  }
  if (!all(lengths == sum)){
    stop("Dimensions of ",.xname," and 'at' must be the same",
         call. = FALSE)
  }
}

.check_verbose <- function(verbose){
  if (!assertive::is_a_bool(verbose)){
    stop("'verbose' must be TRUE or FALSE",
         call. = FALSE)
  }
}

#' @name replaceLetterAt
#' 
#' @title Replacing letters in a nucleotide sequence (or set of nucleotide 
#' sequences) at some specified locations containing nucleotide modifications 
#' 
#' @description 
#' \code{replaceLetterAt} replaces a letter in a \code{\link{ModString}} objects
#' with a new letter. In contrast to \code{\link{modifyNucleotides}} it does not
#' check the letter to be replaced for its identity, it just replaces it and 
#' behaves exactly like the 
#' 
#' @param x a \code{\link{ModString}} or \code{\link{ModStringSet}} object
#' @param at the location where the replacement should be made.
#' 
#' The same input as in \code{\link[Biostrings]{replaceLetterAt}} are expected:
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
#' @param letter The new letters.
#' 
#' The same input as in \code{\link[Biostrings]{replaceLetterAt}} are expected:
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
#' (all(width(letter) == rowSums(at))).
#' @param verbose See \code{\link[Biostrings]{replaceLetterAt}}.
#' 
#' @return the input \code{\link{ModString}} or \code{\link{ModStringSet}}
#' object with the changes applied
#' 
#' @export
#' 
#' @examples
#' # Replacing the last two letters in a ModDNAString
#' seq1 <- ModDNAString("AGTC")
#' seq2 <- replaceLetterAt(seq1,c(3,4),"CT")
#' # Now containg and m3C
#' seq2 <- replaceLetterAt(seq1,c(3,4),ModDNAString("/T"))
#' # Replacing the last two letters in a set of sequences
#' set1 <- ModDNAStringSet(c("AGTC","AGTC"))
#' set2 <- replaceLetterAt(set1,
#'                           matrix(rep(c(FALSE,FALSE,TRUE,TRUE),2),
#'                                  nrow = 2,
#'                                  byrow = TRUE),
#'                           c("CT","CT"))
NULL

#' @rdname replaceLetterAt
#' @export
setMethod(
  "replaceLetterAt",
  signature = "ModString",
  definition = function(x,
                        at,
                        letter,
                        verbose = FALSE){
    .check_verbose(verbose)
    at <- .check_replace_pos_ModString(x,at)
    if (is(letter, "ModString")){
      letter <- as.character(letter)
    }
    else if (!is.character(letter)){
      stop("'letter' must be a ModString object or a character vector",
           call. = FALSE)
    }
    letter <- .convert_letters_to_one_byte_codes(letter,
                                                 modscodec(seqtype(x)))
    .Call2("XString_replace_letter_at",
           x,
           at,
           letter,
           NULL,
           if.not.extending = "replace",
           verbose,
           PACKAGE = "Biostrings")
  }
)

#' @rdname replaceLetterAt
#' @export
setMethod(
  "replaceLetterAt",
  signature = "ModStringSet",
  definition = function(x,
                        at,
                        letter,
                        verbose = FALSE){
    .check_verbose(verbose)
    if (length(x) == 0L){
      stop("'x' has no element")
    }
    .check_replace_pos_ModStringSet(x,
                                    at)
    if (is(letter, "ModStringSet")){
      letter <- as.character(letter)
    }
    else if (!is.character(letter)){
      stop("'letter' must be a ModStringSet object or a character vector")
    }
    .check_letter_ModStringSet(x,at,letter)
    unlisted_x <- unlist(x, use.names = FALSE)
    if(is.list(at)){
      at <- unlist(at)
    } else {
      at <- as.vector(t(at))
    }
    unlisted_ans <- replaceLetterAt(unlisted_x,
                                    at,
                                    letter,
                                    if.not.extending = "replace",
                                    verbose = verbose)
    relist(unlisted_ans, x)
  }
)
