#' @include Modstrings.R
#' @include Modstrings-ModStringCodec.R
NULL

#' @name alphabet
#' 
#' @title Base information for sequence characters of nucleotide strings 
#' containing modifications
#' 
#' @description 
#' The \code{alphabet()}, \code{shortName()} \code{fullName()} and
#' \code{nomenclature()} functions return the letters, names and associated
#' abbreviations for the type of ModString. \code{alphabet()} returns the normal
#' letters and modification letters, whereas \code{shortName()},
#' \code{fullName()} and \code{nomenclature()} return results for modifications
#' only.
#' 
#' @param x a character vector
#' @param baseOnly \code{TRUE} or \code{FALSE} (default): Should the result omit
#' occurances of the letters \code{N.-+}?
#' 
#' @return a character vector.
#' 
#' @examples 
#' alphabet(ModDNAString())
#' shortName(ModDNAString())
#' nomenclature(ModDNAString())
NULL

#' @rdname alphabet
#' @export
setMethod("alphabet", "ModString",
          function(x, baseOnly = FALSE){
            names(Modstrings:::modscodes(seqtype(x),
                                         lettersOnly = TRUE,
                                         baseOnly))
          })
#' @rdname alphabet
#' @export
setMethod("alphabet", "ModStringSet",
          function(x, baseOnly = FALSE){
            names(Modstrings:::modscodes(seqtype(x),
                                         lettersOnly = TRUE,
                                         baseOnly))
          })
#' @rdname alphabet
#' @export
setMethod("shortName", "ModString",
          function(x){
            switch(seqtype(x),
                   ModDNA = MOD_DNA_STRING_CODEC@additionalInfo$short_name,
                   ModRNA = MOD_RNA_STRING_CODEC@additionalInfo$short_name,
                   NULL
            )
          })
#' @rdname alphabet
#' @export
setMethod("shortName", "ModStringSet",
          function(x){
            switch(seqtype(x),
                   ModDNA = MOD_DNA_STRING_CODEC@additionalInfo$short_name,
                   ModRNA = MOD_RNA_STRING_CODEC@additionalInfo$short_name,
                   NULL
            )
          })
#' @rdname alphabet
#' @export
setMethod("fullName", "ModString",
          function(x){
            switch(seqtype(x),
                   ModDNA = MOD_DNA_STRING_CODEC@additionalInfo$name,
                   ModRNA = MOD_RNA_STRING_CODEC@additionalInfo$name,
                   NULL
            )
          })
#' @rdname alphabet
#' @export
setMethod("fullName", "ModStringSet",
          function(x){
            switch(seqtype(x),
                   ModDNA = MOD_DNA_STRING_CODEC@additionalInfo$name,
                   ModRNA = MOD_RNA_STRING_CODEC@additionalInfo$name,
                   NULL
            )
          })
#' @rdname alphabet
#' @export
setMethod("nomenclature", "ModString",
          function(x){
            switch(seqtype(x),
                   ModDNA = MOD_DNA_STRING_CODEC@additionalInfo$nc,
                   ModRNA = MOD_RNA_STRING_CODEC@additionalInfo$nc,
                   NULL
            )
          })
#' @rdname alphabet
#' @export
setMethod("nomenclature", "ModStringSet",
          function(x){
            switch(seqtype(x),
                   ModDNA = MOD_DNA_STRING_CODEC@additionalInfo$nc,
                   ModRNA = MOD_RNA_STRING_CODEC@additionalInfo$nc,
                   NULL
            )
          })


# derived from Biostrings/R/seqtype.R ------------------------------------------


.ModDNAorRNAcodes <- function(mods,
                              base_codes,
                              lettersOnly,
                              baseOnly,
                              col_name = c("abbrev",
                                           "short_name",
                                           "nc")) {
  # remove empty letters. this is four neutrality against currently 
  # unsupported modifications. 
  # for more details have a look add the ModStringCodec class
  mods <- mods[mods$abbrev != "",]
  #
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

base_class_name <- function(x) paste(seqtype(x), "String", sep="")

modscodec <- function(x){
  switch(x,
         ModDNA = MOD_DNA_STRING_CODEC,
         ModRNA = MOD_RNA_STRING_CODEC,
         NULL
  )
}

modscodes <- function(x,
                      lettersOnly = FALSE,
                      baseOnly = FALSE){
  switch(x,
         ModDNA = .ModDNAorRNAcodes(MOD_DNA_BASE_CODES,
                                    Biostrings:::DNA_BASE_CODES,
                                    lettersOnly,
                                    baseOnly),
         ModRNA = .ModDNAorRNAcodes(MOD_RNA_BASE_CODES,
                                    Biostrings:::RNA_BASE_CODES,
                                    lettersOnly,
                                    baseOnly),
         0:255
  )
}

modsshortnames <- function(x){
  switch(x,
         ModDNA = .ModDNAorRNAcodes(MOD_DNA_BASE_CODES,
                                    c(),
                                    TRUE,
                                    TRUE,
                                    "short_name"),
         ModRNA = .ModDNAorRNAcodes(MOD_RNA_BASE_CODES,
                                    c(),
                                    TRUE,
                                    TRUE,
                                    "short_name"),
         NULL
  )
}

modsnomenclature <- function(x){
  switch(x,
         ModDNA = .ModDNAorRNAcodes(MOD_DNA_BASE_CODES,
                                    c(),
                                    TRUE,
                                    TRUE,
                                    "nc"),
         ModRNA = .ModDNAorRNAcodes(MOD_RNA_BASE_CODES,
                                    c(),
                                    TRUE,
                                    TRUE,
                                    "nc"),
         NULL
  )
}

comparable_seqtypes <- function(seqtype1,
                                seqtype2){
  is_nucleo1 <- seqtype1 %in% c("ModDNA","ModRNA","DNA","RNA")
  is_nucleo2 <- seqtype2 %in% c("ModDNA","ModRNA","DNA","RNA")
  is_same_type <- (seqtype1 %in% c("ModDNA","DNA") & 
                     seqtype2 %in% c("ModDNA","DNA")) |
    (seqtype1 %in% c("ModRNA","RNA") & 
       seqtype2 %in% c("ModRNA","RNA"))
  is_nucleo1 == is_nucleo2 & is_same_type
}
