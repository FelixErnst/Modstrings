#' @include Modstrings.R
NULL

additional_base_codes <- c(N = 15L, `-` = 16L, `+` = 32L, `.` = 64L)

.load_mod_dictionary <- function(file){
  con <- file(file)
  file <- stringr::str_split(readLines(con,
                                       encoding = "UTF-8"),
                             "\t")
  close(con)
  # remove commented lines
  file <- file[!substr(lapply(file,"[",1),1,1) %in% c("#","")]
  #
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

MOD_DNA_BASE_CODES <- .load_mod_dictionary("inst/extdata/Mod_DNA_codes.txt")
MOD_RNA_BASE_CODES <- .load_mod_dictionary("inst/extdata/Mod_RNA_codes.txt")

################################################################################
# ModStringCodec is used to convert a strubg with multi byte letters into string
# with one byte letters. The resulting object can then be used like a BString.
# 
# However, this requires the conversion during printing or export. --> see below
# for conversion functions. These function are used by the Mod*String class
################################################################################

setClass("ModStringCodec",
         slots = c(letters = "character",
                   oneByteCodes = "character",
                   conversion = "logical",
                   originatingBase = "character",
                   values = "integer",
                   lettersEscaped = "character",
                   oneByteCodesEscaped = "character",
                   lettersNeedEscape = "logical",
                   oneByteCodesNeedEscape = "logical",
                   additionalInfo = "DataFrame"))

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
.str_replace_all_fixed_custom <- function(string,
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
.str_replace_all_regex_custom <- function(string,
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

.new_ModStringCodec <- function(base_codes, biostrings_base_codes){
  letters <- base_codes$abbrev
  oneByteCodes <- base_codes$oneByteLetter
  orig_base <- base_codes$orig_base
  values <- base_codes$value
  extra_letters <- biostrings_base_codes
  additionalInfo <- base_codes[,c("name","short_name","nc","orig_base","abbrev")]
    lengths <- unique(c(length(letters),
                        length(oneByteCodes),
                        length(values),
                        length(orig_base)))
  if(length(lengths) != 1){
    stop("ModStringCodec: Input do not have the same length.")
  }
  # remove empty letters. this is four neutrality against currently 
  # unsupported modifications. However they can be part of the
  # additionalInfo, which is used for the construction of the 
  # sanitization dictionaries
  f <- letters == ""
  letters <- letters[!f]
  oneByteCodes <- oneByteCodes[!f]
  values <- values[!f]
  orig_base <- orig_base[!f]
  # 
  letters <- c(letters,names(extra_letters))
  oneByteCodes <- c(oneByteCodes,names(extra_letters))
  originatingBase <- c(orig_base,names(extra_letters))
  values <- c(values,unname(extra_letters))
  # originating base must be in the extra_letter or empty
  if(!all(originatingBase %in% c(names(extra_letters),""))){
    stop("Not all originating bases are in the extra letters.")
  }
  # order based on the values in ascending order
  f <- order(values)
  letters <- letters[f]
  oneByteCodes <- oneByteCodes[f]
  originatingBase <- originatingBase[f]
  values <- values[f]
  # escape necessary values:
  lettersEscaped <- .escape_special_charactes(letters)
  lettersNeedEscape <- letters != lettersEscaped
  oneByteCodesEscaped <- .escape_special_charactes(oneByteCodes)
  oneByteCodesNeedEscape <- oneByteCodes != oneByteCodesEscaped
  # check which letters need conversion
  # and control input
  conversion <- vapply(letters,
                       function(l){
                         length(charToRaw(l)) > 1
                       },
                       logical(1))
  checkConversionValid <- vapply(c(letters[!which(conversion)],
                                   oneByteCodes[which(conversion)]),
                                 function(l){
                                   length(charToRaw(l)) > 1
                                 },
                                 logical(1))
  if(any(checkConversionValid)){
    stop("Not all letters have a proper one byte conversion.")
  }
  new("ModStringCodec", letters = letters, oneByteCodes = oneByteCodes, 
      conversion = conversion, originatingBase = originatingBase,
      values = values, lettersEscaped = lettersEscaped,
      oneByteCodesEscaped = oneByteCodesEscaped, 
      lettersNeedEscape = lettersNeedEscape,
      oneByteCodesNeedEscape = oneByteCodesNeedEscape, 
      additionalInfo = additionalInfo)
}

MOD_DNA_STRING_CODEC <- .new_ModStringCodec(MOD_DNA_BASE_CODES,
                                            c(Biostrings:::DNA_BASE_CODES,
                                              additional_base_codes))
MOD_RNA_STRING_CODEC <- .new_ModStringCodec(MOD_RNA_BASE_CODES,
                                            c(Biostrings:::RNA_BASE_CODES,
                                              additional_base_codes))
