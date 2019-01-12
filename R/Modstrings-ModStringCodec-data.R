#' @include Modstrings.R

################################################################################
# This file contains helper functions only, which are used to generate and 
# maintain the dictionaries
################################################################################

# loading/writing tab delimeted file

# .write_to_tab <- function(df,file,header = FALSE){
#   if(header) {
#     input <- paste0(colnames(df),collapse = "\t")
#   } else {
#     input <- c()
#   }
#   input <- c(input,
#              unlist(lapply(seq_len(nrow(df)),
#                            function(i){
#                              paste0(df[i,,drop=TRUE],collapse = "\t")
#                            })))
#   writeLines(input,
#              con = file,
#              useBytes = TRUE)
# }
# .load_from_tab_to_DataFrame <- function(file){
#   as(read.delim(file,
#                 encoding = "UTF-8",
#                 sep = "\t",
#                 quote = "",
#                 stringsAsFactors = FALSE),"DataFrame")
# }

# Dictionaries
# for tRNAdb sanitization ------------------------------------------------------

# MOD_RNA_DICT_TRNADB <- read.delim("inst/extdata/tRNAdb_modifications.txt",
#                        header = FALSE,
#                        quote = "",
#                        stringsAsFactors = FALSE)
# colnames(MOD_RNA_DICT_TRNADB) <- c("rnamods_abbrev","short_name","name")

# for modomics sanitization ----------------------------------------------------

# MOD_RNA_DICT_MODOMICS <- read.delim("inst/extdata/modomics_modifications.txt",
#                                     header = FALSE,
#                                     quote = "",
#                                     stringsAsFactors = FALSE)
# colnames(MOD_RNA_DICT_TRNADB) <- c("rnamods_abbrev","short_name","name")

################################################################################

# Working on the alphabets

# ModDNA alphabet -----------------------------------------------------------------

# modsDNA <- as(read.delim("inst/extdata/DNAmod_modifications.txt",
#                       encoding = "UTF-8", sep = " ", quote = "",
#                       stringsAsFactors = FALSE),"DataFrame")
# modsDNA$nc <- modsDNA$short_name

# ModRNA alphabet -----------------------------------------------------------------

# modsDNA <- as(read.delim("inst/extdata/DNAmod_modifications.txt",
#                       encoding = "UTF-8", sep = " ", quote = "",
#                       stringsAsFactors = FALSE),"DataFrame")
# modsDNA$nc <- modsDNA$short_name

# fixes for certain characters, which are not compatible

# mods[mods$short_name == "xX","rnamods_abbrev"] <- "÷"
# mods[mods$short_name == "xU","rnamods_abbrev"] <- "Ü"
# mods[mods$short_name == "mcmo5Um","rnamods_abbrev"] <- "Ä"
# mods[mods$short_name == "m5Um","rnamods_abbrev"] <- "¤"
# mods[mods$short_name == "f5C","rnamods_abbrev"] <- "×"
# mods[mods$short_name == "N","rnamods_abbrev"] <- "" 

# generating value and oneByteLetters in a sensible way for avoiding clashes
# and other things
# .add_oneByteLetters_to_moficiation_codes <- function(mods,
#                                                      base_codes){
#   mods$value <- 0L
#   mods$oneByteLetter <- 0L
#   # avoid clash with predefined codes
#   # assumes that non are part of mods
#   codes <- c(base_codes,additional_base_codes)
#   browser()
#   if(any(names(codes) %in% mods$abbrev)){
#     stop("Abbreviations cannot be the following characters: '",
#          paste(names(codes), collapse = "','"),
#          "'")
#   }
#   n <- nrow(mods) + length(codes)
#   # if less codes are needed than the highest of the base codes,
#   # generate less codes
#   n <- n - length(codes[codes > n])
#   ints <- as.integer(seq(n))
#   ints <- ints[!(ints %in% codes)]
#   # add integer value to codes
#   mods[mods$value == 0L,]$value <- ints
#   # since modifications can have more than one byte per character, create a
#   # dynamic letter to be used internally. This requires additional conversion
#   # for printing and saving, but not for computing
#   ###########
#   # 1. set the oneByteLetter from the abbreviations and check which are to
#   # long
#   mods$oneByteLetter <- mods$abbrev
#   f <- vapply(mods$oneByteLetter,
#               function(obl){
#                 length(charToRaw(obl)) > 1
#               },
#               logical(1))
#   f_neg <- which(!f) # not too long
#   f <- which(f) # too long
#   # 2. get all avialble one byte letters
#   oneByteLetter <- unlist(lapply(seq_len(255),
#                                  function(i){rawToChar(as.raw(i))}))
#   # 3. remove the default codes and the ones already present
#   oneByteLetter <- oneByteLetter[
#     !duplicated(oneByteLetter) &
#       !(oneByteLetter %in% names(codes)) &
#       !(oneByteLetter %in% mods$abbrev[f_neg])]
#   # 4. if by any chance not enough one byte letters are present
#   # annoy me immediately
#   if(length(oneByteLetter) < nrow(mods)){
#     stop("Somebody broke it. Cannot handle that many different ",
#          "modifications.",
#          call. = FALSE)
#   }
#   # 5. exchange the letters, which are to long, from the available pool of one
#   # byte letters
#   mods$oneByteLetter[f] <- oneByteLetter[seq_len(length(f))]
#   mods
# }

# convert oneByteLetter back to oneByteInteger

# .convert_oneByteLetters_to_oneByteInteger <- function(mods){
#   mods$oneByteInteger <- 
#     as.integer(unlist(lapply(mods$oneByteLetter,charToRaw)))
#   mods$oneByteLetter <- NULL
#   mods
# }

# 
# MOD_DNA_BASE_CODES <- .add_oneByteLetters_to_moficiation_codes(
#   DNAmods,
#   Biostrings:::DNA_BASE_CODES)
# MOD_RNA_BASE_CODES <- .add_oneByteLetters_to_moficiation_codes(
#   RNAmods,
#   Biostrings:::RNA_BASE_CODES)
# MOD_DNA_BASE_CODES <- 
#   .convert_oneByteLetters_to_oneByteInteger(MOD_DNA_BASE_CODES)
# MOD_RNA_BASE_CODES <- 
#   .convert_oneByteLetters_to_oneByteInteger(MOD_RNA_BASE_CODES)
