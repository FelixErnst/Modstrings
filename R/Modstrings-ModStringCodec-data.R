#' @include Modstrings.R

################################################################################
# This file contains helper functions only, which are used to generate and 
# maintain the dictionaries
################################################################################

# MOD_RNA_DICT_TRNADB <- read.delim("inst/extdata/tRNAdb_modifications.txt",
#                        header = FALSE,
#                        quote = "Ü",
#                        stringsAsFactors = FALSE)
# colnames(MOD_RNA_DICT_TRNADB) <- c("rnamods_abbrev","short_name","name")

################################################################################
# preload modification data and remove indistinguishable modifications
# This data is UTF-8 encoded. Every change to the modification abbreviations
# has to UTF-8 encoded as well
# modsDNA <- as(read.delim("inst/extdata/DNAmod_modifications.txt",
#                       encoding = "UTF-8", sep = " ", quote = "",
#                       stringsAsFactors = FALSE),"DataFrame")
# modsDNA$nc <- modsDNA$short_name
# load("data/mods.rda")
# load("data/modsDNA.rda")
# #
# # # fix some abbreviations
# mods[mods$short_name == "xX","rnamods_abbrev"] <- "÷" # @ to ¶ for fastq
# mods[mods$short_name == "xU","rnamods_abbrev"] <- "Ü"
# mods[mods$short_name == "mcmo5Um","rnamods_abbrev"] <- "Ä"
# mods[mods$short_name == "m5Um","rnamods_abbrev"] <- "¤"
# mods[mods$short_name == "f5C","rnamods_abbrev"] <- "×"
# mods[mods$short_name == "N","rnamods_abbrev"] <- "" # available from base
# 
# mods <- mods[!(mods$rnamods_abbrev %in%
#               unique(mods$rnamods_abbrev[duplicated(mods$rnamods_abbrev)])),]
# modomicColNames <- c("name",
#                      "short_name",
#                      "new_nomenclature",
#                      "originating_base",
#                      "rnamods_abbrev",
#                      "html_abbrev")
# modColNames <- c("name",
#                  "short_name",
#                  "nc",
#                  "orig_base",
#                  "abbrev",
#                  "html_abbrev")
# mods[mods$originating_base == "preQ0base","originating_base"] <- "G"
# mods <- mods[!stringi::stri_detect_fixed(mods[,"new_nomenclature"],"(base)"),]
# ##############################################################################
# # split into RNA and DNA modifications
# RNAmods <- mods[!(mods$originating_base  %in% c("pppN")) &
#                   !(mods$rnamods_abbrev  %in% c("","none",".")) &
#                   # !grepl("unknown",mods$name) &
#                   !(mods$name %in% c("adenosine",
#                                      "uridine",
#                                      "cytidine",
#                                      "guanosine")),
#                 modomicColNames]
# DNAmods <- modsDNA
# colnames(DNAmods) <- modColNames[c(1,2,5,4,3)]
# colnames(RNAmods) <- modColNames
# 
# ##############################################################################
# .add_oneByteLetters_to_moficiation_codes <- function(mods,
#                                                      base_codes){
#   mods$value <- 0L
#   mods$oneByteLetter
#   # avoid clash with predefined codes
#   # assumes that non are part of mods
#   codes <- c(base_codes,additional_base_codes)
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
#   f_neg <- which(!f)
#   f <- which(f)
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
# 
# MOD_DNA_BASE_CODES <- .add_oneByteLetters_to_moficiation_codes(
#   DNAmods,
#   Biostrings:::DNA_BASE_CODES)
# MOD_RNA_BASE_CODES <- .add_oneByteLetters_to_moficiation_codes(
#   RNAmods,
#   Biostrings:::RNA_BASE_CODES)
# usethis::use_data(MOD_DNA_BASE_CODES, overwrite = TRUE)
# usethis::use_data(MOD_RNA_BASE_CODES, overwrite = TRUE)
