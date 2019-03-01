#' @include Modstrings.R
NULL

#' @name sanitizeInput
#' @aliases sanitizeFromModomics sanitizeFromtRNAdb
#'
#' @title Sanitize input strings for use with ModString classes
#' 
#' @description 
#' Since the one letter nomenclature for RNA and DNA modification differs
#' depending on the source, a translation to a common alphabet is necessary.
#' 
#' \code{sanitizeInput} exchanges based on a dictionary. The dictionary is
#' expected to be a \code{DataFrame} with two columns, \code{mods_abbrev} and
#' \code{short_name}. Based on the \code{short_name} the characters from in the 
#' input are converted from values of \code{mods_abbrev} into the the ones
#' from \code{alphabet}.
#' 
#' Only different values will be searched for and exchanged.
#' 
#' \code{sanitizeFromModomics} and \code{sanitizeFromtRNAdb} use a predefined
#' dictionary, which is builtin.
#' 
#' @param input a \code{character} vector, which should be converted
#' @param dictionary a DataFrame containing at least two columns 
#' \code{mods_abbrev} and \code{short_name}. From this a dictionary table is
#' contructed for exchaning old to new letters.
#'
#' @return the modified \code{character} vector compatible for constructing a
#' \code{ModString} object.
#' @export
#'
#' @examples
#' # Modomics
#' chr <- "AGC@"
#' # Error since the @ is not in the alphabet
#' \dontrun{
#' seq <- ModRNAString(chr)
#' }
#' seq <- ModRNAString(sanitizeFromModomics(chr))
#' 
#' # tRNAdb
#' chr <- "AGC+"
#' # No error but the + has a different meaning in the alphabet
#' \dontrun{
#' seq <- ModRNAString(chr)
#' }
#' seq <- ModRNAString(sanitizeFromtRNAdb(chr))
NULL

load("data/MOD_RNA_DICT_MODOMICS.rda")
load("data/MOD_RNA_DICT_TRNADB.rda")

# this is required since toupper() screws up the encoding of special characters
.add_lower_to_upper_case <- function(dict)
{
  caseDict <- DataFrame(old = c("a","g","c","t","u"),
                    sn = c("A","G","C","T","U"),
                    new = c("A","G","C","T","U"))
  caseDict <- caseDict[!(caseDict$old %in% dict$old),]
  dict <- rbind(dict,
                caseDict)
  dict
}

.construct_translation_table <- function(dict, seqtype)
{
  if(!is(dict,"DataFrame")){
    stop("Dictionary must be a DataFrame with at least two columns ",
         "'mods_abbrev' and 'short_name'. 'short_name' must match the ",
         "nomenclature returned by nomenclature().",
         call. = FALSE)
  }
  if("rnamods_abbrev" %in% colnames(dict)){
    colnames <- colnames(dict)
    colnames[colnames == "rnamods_abbrev"] <- "mods_abbrev"
    colnames(dict) <- colnames
  }
  if(seqtype == "ModRNA"){
    internalDict <- additionalInfo(MOD_RNA_STRING_CODEC)
  } else {
    internalDict <- additionalInfo(MOD_DNA_STRING_CODEC)
  }
  dict <- dict[dict$short_name %in% internalDict$short_name,c("mods_abbrev",
                                                              "short_name")]
  colnames(dict) <- c("old","sn")
  dict$new <- internalDict[match(dict$sn,internalDict$short_name),"abbrev"]
  dict$old <- as.character(dict$old)
  dict$sn <- as.character(dict$sn)
  dict$new <- as.character(dict$new)
  dict <- dict[dict$old != dict$new,]
  dict <- .add_lower_to_upper_case(dict)
  dict
}

#' @rdname sanitizeInput
#' @export
sanitizeInput <- function(input, dictionary)
{
  if(is.list(input)){
    input <- unlist(input)
  }
  if(!is.character(input)){
    stop("Input has to be a 'character' vector or a list containing only ",
         "'character' vectors.",
         call. = FALSE)
  }
  dict <- .construct_translation_table(dictionary,"ModRNA")
  if(nrow(dict) > 0L){
    names <- names(input)
    if(stringi::stri_enc_get() == "UTF-8"){
      input <- vapply(input,
                      .str_replace_all_fixed_custom,
                      character(1),
                      dict$old,
                      dict$new,
                      USE.NAMES = FALSE)
    } else {
      input <- vapply(input,
                      .str_replace_all_regex_custom,
                      character(1),
                      .escape_special_charactes(dict$old),
                      dict$new,
                      USE.NAMES = FALSE)
    }
    names(input) <- names
  }
  input
}

#' @rdname sanitizeInput
#' @export
sanitizeFromModomics <- function(input){
  sanitizeInput(input, MOD_RNA_DICT_MODOMICS)
}

#' @rdname sanitizeInput
#' @export
sanitizeFromtRNAdb <- function(input){
  sanitizeInput(input, MOD_RNA_DICT_TRNADB)
}
