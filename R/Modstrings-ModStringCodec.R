#' @include Modstrings.R
NULL

additional_base_codes <- c(N = 15L, `-` = 16L, `+` = 32L, `.` = 64L)

load_mod_dictionary <- function(file){
  con <- file(file)
  file <- stringr::str_split(readLines(con,
                                       encoding = "UTF-8"),
                             "\t")
  close(con)
  m <- t(unname(data.frame(file,
                           stringsAsFactors = FALSE)))
  colnames(m) <- m[1,]
  m <- m[-1,]
  df <- DataFrame(m)
  df$name <- as.character(df$name)
  df$short_name <- as.character(df$short_name)
  df$abbrev <- as.character(df$abbrev)
  df$orig_base <- as.character(df$orig_base)
  df$nc <- as.character(df$nc)
  df$value <- as.integer(m[,"value"])
  df$oneByteInteger <- as.integer(m[,"oneByteInteger"])
  df$oneByteLetter <- unlist(lapply(lapply(as.integer(df$oneByteInteger),
                                           as.raw),
                                    rawToChar))
  df[,!(colnames(df) %in% c("oneByteInteger"))]
}

#load the dictionaries from file

MOD_DNA_BASE_CODES <- load_mod_dictionary("inst/extdata/Mod_DNA_codes.txt")
MOD_RNA_BASE_CODES <- load_mod_dictionary("inst/extdata/Mod_RNA_codes.txt")

################################################################################
# ModStringCodec is used to convert the ModString with multi byte letters
# into string with one byte letters
# This can then be used as BString. However, this requires the conversion
# during printing or export. --> see below for conversion functions
################################################################################

setClass("ModStringCodec",
         slots = c(
           letters = "character",
           oneByteCodes = "character",
           conversion = "logical",
           originatingBase = "character",
           values = "integer",
           lettersEscaped = "character",
           oneByteCodesEscaped = "character",
           lettersNeedEscape = "logical",
           oneByteCodesNeedEscape = "logical",
           additionalInfo = "DataFrame"
         ))
setMethod("initialize",
          "ModStringCodec",
          function(.Object,
                   letters,
                   oneByteCodes,
                   orig_base,
                   values,
                   extra_letters,
                   additionalInfo){
            lengths <- unique(c(length(letters),
                                length(oneByteCodes),
                                length(values),
                                length(orig_base)))
            if(length(lengths) != 1){
              stop("ModStringCodec: Input do not have the same length.")
            }
            .Object@letters <- c(letters,names(extra_letters))
            .Object@oneByteCodes <- c(oneByteCodes,names(extra_letters))
            .Object@originatingBase <- c(orig_base,names(extra_letters))
            .Object@values <- c(values,unname(extra_letters))
            .Object@additionalInfo <- additionalInfo
            # originating base must be in the extra_letter or empty
            if(!all(.Object@originatingBase %in% c(names(extra_letters),""))){
              stop("Not all originating bases are in the extra letters.")
            }
            # order
            .Object@letters <- .Object@letters[order(.Object@values)]
            .Object@oneByteCodes <- .Object@oneByteCodes[order(.Object@values)]
            .Object@originatingBase <- 
              .Object@originatingBase[order(.Object@values)]
            .Object@values <- .Object@values[order(.Object@values)]
            # escape necessary values:
            # -
            # -
            .Object@lettersEscaped <- 
              .escape_special_charactes(.Object@letters)
            .Object@lettersNeedEscape <- 
              .Object@letters != .Object@lettersEscaped
            #
            .Object@oneByteCodesEscaped <- 
              .escape_special_charactes(.Object@oneByteCodes)
            .Object@oneByteCodesNeedEscape <- 
              .Object@oneByteCodes != .Object@oneByteCodesEscaped
            # check which letters need conversion
            # and control input
            f <- vapply(.Object@letters,
                        function(l){
                          length(charToRaw(l)) > 1
                        },
                        logical(1))
            .Object@conversion <- f
            f <- vapply(c(.Object@letters[!which(f)],
                          .Object@oneByteCodes[which(f)]),
                        function(l){
                          length(charToRaw(l)) > 1
                        },
                        logical(1))
            if(any(f)){
              stop("Not all letters have a proper one byte conversion.")
            }
            .Object
          }
)

#' @importFrom stringr str_detect
.escape_special_charactes <- function(x){
  special_chars <- paste0("\\~\\!\\#\\$\\%\\^\\&\\*\\{\\}\\+\\:",
                          "\\\\\"\\?\\,\\.\\/\\'\\[\\]\\-\\|")
  f <- which(stringr::str_detect(x,paste0("[:punct:]|[",special_chars,"]")))
  x[f] <- paste0("\\",x[f])
  x
}

################################################################################
# additional sequence conversion functions, since some characters are longer
# than one byte. They need to encoded first using the information in the 
# modification base code DataFrame.
####
# The incoming string does need to be UTF-8 encoded and the outgoing must be
# ASCII one byte letters
################################################################################

#' @importFrom stringi stri_enc_toutf8 stri_enc_toascii stri_enc_get 
#' stri_enc_tonative
.convert_letters_to_one_byte_codes <- function(string,
                                               codec){
  string <- stringi::stri_enc_toutf8(string)
  .check_for_invalid_letters(string,
                             codec)
  obc_string <- .str_replace_all_regex_custom(
    string,
    codec@lettersEscaped[codec@conversion],
    codec@oneByteCodes[codec@conversion])
  if(stringi::stri_enc_get() == "UTF-8"){
    return(obc_string)
  }
  return(stringi::stri_enc_tonative(obc_string))
}

.check_for_invalid_letters <- function(string,
                                       codec){
  letters_in_string <- unique(strsplit(string,"")[[1]])
  if(any(!(letters_in_string %in% codec@letters))){
    print(paste(
      letters_in_string[!(letters_in_string %in% codec@letters)],
      collapse = ""))
    stop("Invalid character(s) - see above",
         call. = FALSE)
  }
}

.convert_one_byte_codes_to_letters <- function(obc_string,
                                               codec){
  obc_string <- stringi::stri_enc_toutf8(obc_string)
  if(stringi::stri_enc_get() == "UTF-8"){
    string <- .str_replace_all_fixed_custom(
      obc_string,
      codec@oneByteCodes[codec@conversion],
      codec@letters[codec@conversion])
  } else {
    string <- .str_replace_all_regex_custom(
      obc_string,
      codec@oneByteCodesEscaped[codec@conversion],
      codec@letters[codec@conversion])
  }
  return(string) 
}

#' @importFrom stringi stri_locate_all_fixed stri_sub
.str_replace_all_fixed_custom_c <- function(string,
                                            pattern,
                                            replacement){
  locations <- stringi::stri_locate_all_fixed(string,
                                              pattern,
                                              omit_no_match = TRUE)
  f <- which(!vapply(locations,function(l){nrow(l) == 0},logical(1)))
  # Currently now idea how to avoid the loops
  for(i in f){
    loc <- locations[[i]]
    for(j in seq_len(nrow(loc))){
      stringi::stri_sub(string,loc[j,"start"],loc[j,"end"]) <- 
        replacement[i]
    }
  }
  string
}
#' @importFrom stringi stri_locate_all_regex stri_sub
.str_replace_all_regex_custom_c <- function(string,
                                            pattern,
                                            replacement){
  locations <- stringi::stri_locate_all_regex(string,
                                              pattern,
                                              omit_no_match = TRUE)
  f <- which(!vapply(locations,function(l){nrow(l) == 0},logical(1)))
  # Currently now idea how to avoid the loops
  for(i in f){
    loc <- locations[[i]]
    for(j in seq_len(nrow(loc))){
      stringi::stri_sub(string,loc[j,"start"],loc[j,"end"]) <- 
        replacement[i]
    }
  }
  string
}
.str_replace_all_regex_custom <- 
  compiler::cmpfun(.str_replace_all_regex_custom_c)
.str_replace_all_fixed_custom <- 
  compiler::cmpfun(.str_replace_all_fixed_custom_c)


.convert_letters_to_originating_base <- function(string,
                                                 codec){
  string <- stringi::stri_enc_toutf8(string)
  orig_string <- .str_replace_all_regex_custom(
    string,
    codec@lettersEscaped,
    codec@originatingBase)
  if(stringi::stri_enc_get() == "UTF-8"){
    return(orig_string)
  }
  return(stringi::stri_enc_tonative(orig_string))
}

.convert_one_byte_codes_to_originating_base <- function(obc_string,
                                                        codec){
  obc_string <- stringi::stri_enc_toutf8(obc_string)
  if(stringi::stri_enc_get() == "UTF-8"){
    orig_string <- .str_replace_all_fixed_custom(
      obc_string,
      codec@oneByteCodes,
      codec@originatingBase)
  } else {
    orig_string <- .str_replace_all_regex_custom(
      obc_string,
      codec@oneByteCodesEscaped,
      codec@originatingBase)
  }
  return(orig_string) 
}

################################################################################
# get Biostring like base codes, alphabet and Codec object
# the codec object is not inherited from Biostrings package, but is
# used for one byte conversion only

.ModStringCodec <- function(base_codes,
                            biostrings_base_codes){
  new("ModStringCodec", 
      base_codes$abbrev,
      base_codes$oneByteLetter,
      base_codes$orig_base,
      base_codes$value,
      biostrings_base_codes,
      base_codes[,c("name","short_name","nc","orig_base","abbrev")])
}

.ModDNAorRNAcodes <- function(mods,
                              base_codes,
                              lettersOnly,
                              baseOnly,
                              col_name = c("abbrev",
                                           "short_name",
                                           "nc")) {
  col_name <- match.arg(col_name)
  codes <- mods$value
  if(is.null(mods$oneByteLetter) || lettersOnly){
    names(codes) <- mods[,col_name]
  } else {
    names(codes) <- mods$oneByteLetter
  }
  if(baseOnly){
    return(c(base_codes,codes))
  } else {
    return(c(base_codes,additional_base_codes,codes))
  }
}
MOD_DNA_STRING_CODEC <- .ModStringCodec(MOD_DNA_BASE_CODES,
                                        c(Biostrings:::DNA_BASE_CODES,
                                          additional_base_codes))
MOD_RNA_STRING_CODEC <- .ModStringCodec(MOD_RNA_BASE_CODES,
                                        c(Biostrings:::RNA_BASE_CODES,
                                          additional_base_codes))

MOD_DNA_ALPHABET <- names(.ModDNAorRNAcodes(MOD_DNA_BASE_CODES,
                                            Biostrings:::DNA_BASE_CODES,
                                            TRUE,
                                            FALSE))
MOD_RNA_ALPHABET <- names(.ModDNAorRNAcodes(MOD_RNA_BASE_CODES,
                                            Biostrings:::RNA_BASE_CODES,
                                            TRUE,
                                            FALSE))