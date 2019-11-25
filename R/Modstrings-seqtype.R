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
#' @param x a \code{ModString} or \code{ModStringSet} object
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
            names(xscodes(x, baseOnly = baseOnly, multiByteLetterNames = TRUE))
          })
#' @rdname alphabet
#' @export
setMethod("alphabet", "ModStringSet",
          function(x, baseOnly = FALSE){
            names(xscodes(x, baseOnly = baseOnly, multiByteLetterNames = TRUE))
          })
#' @rdname alphabet
#' @export
setMethod("shortName", "ModString",
          function(x){
            switch(seqtype(x),
                   ModDNA = additionalInfo(MOD_DNA_STRING_CODEC)$short_name,
                   ModRNA = additionalInfo(MOD_RNA_STRING_CODEC)$short_name,
                   NULL
            )
          })
#' @rdname alphabet
#' @export
setMethod("shortName", "ModStringSet",
          function(x){
            switch(seqtype(x),
                   ModDNA = additionalInfo(MOD_DNA_STRING_CODEC)$short_name,
                   ModRNA = additionalInfo(MOD_RNA_STRING_CODEC)$short_name,
                   NULL
            )
          })
#' @rdname alphabet
#' @export
setMethod("fullName", "ModString",
          function(x){
            switch(seqtype(x),
                   ModDNA = additionalInfo(MOD_DNA_STRING_CODEC)$name,
                   ModRNA = additionalInfo(MOD_RNA_STRING_CODEC)$name,
                   NULL
            )
          })
#' @rdname alphabet
#' @export
setMethod("fullName", "ModStringSet",
          function(x){
            switch(seqtype(x),
                   ModDNA = additionalInfo(MOD_DNA_STRING_CODEC)$name,
                   ModRNA = additionalInfo(MOD_RNA_STRING_CODEC)$name,
                   NULL
            )
          })
#' @rdname alphabet
#' @export
setMethod("nomenclature", "ModString",
          function(x){
            switch(seqtype(x),
                   ModDNA = additionalInfo(MOD_DNA_STRING_CODEC)$nc,
                   ModRNA = additionalInfo(MOD_RNA_STRING_CODEC)$nc,
                   NULL
            )
          })
#' @rdname alphabet
#' @export
setMethod("nomenclature", "ModStringSet",
          function(x){
            switch(seqtype(x),
                   ModDNA = additionalInfo(MOD_DNA_STRING_CODEC)$nc,
                   ModRNA = additionalInfo(MOD_RNA_STRING_CODEC)$nc,
                   NULL
            )
          })


# derived from Biostrings/R/seqtype.R ------------------------------------------


.ModDNAorRNAcodes <- function(mods, base_codes, baseOnly = FALSE,
                              multiByteLetterNames = FALSE, 
                              oneByteIntegerValue = TRUE,
                              col_name = c("abbrev", "short_name", "nc"))
{
  # remove empty letters. this is four neutrality against currently 
  # unsupported modifications. 
  # for more details have a look add the ModStringCodec class
  mods <- mods[mods$abbrev != "",]
  #
  col_name <- match.arg(col_name)
  codes <- mods$value
  if(multiByteLetterNames){
    names(codes) <- mods[,col_name]
  } else {
    names(codes) <- mods$oneByteLetter
  }
  if(baseOnly){
    ans <- c(base_codes,codes)
  } else {
    ans <- c(base_codes,additional_base_codes,codes)
  }
  if(oneByteIntegerValue){
    if(baseOnly){
      values <- c(names(base_codes),mods$oneByteLetter)
    } else {
      values <- c(names(c(base_codes,additional_base_codes)),mods$oneByteLetter)
    }
    ans[] <- as.integer(unlist(lapply(values, charToRaw)))
  }
  ans
}

setClassUnion("Modstrings", c("ModString","ModStringSet"))

#' @importFrom Biostrings xscodes
setMethod("xscodes","Modstrings",
          function(x, baseOnly = FALSE, multiByteLetterNames = FALSE){
            if (!assertive::is_a_bool(baseOnly)){
              stop("'baseOnly' must be TRUE or FALSE")
            }
            if (!assertive::is_a_bool(multiByteLetterNames)){
              stop("'baseOnly' must be TRUE or FALSE")
            }
            switch(seqtype(x),
                   ModDNA = .ModDNAorRNAcodes(MOD_DNA_BASE_CODES,
                                              .DNA_BASE_CODES,
                                              baseOnly,
                                              multiByteLetterNames),
                   ModRNA = .ModDNAorRNAcodes(MOD_RNA_BASE_CODES,
                                              .RNA_BASE_CODES,
                                              baseOnly,
                                              multiByteLetterNames),
                   0:255
            )
          })

modscodec <- function(x)
{
  switch(x,
         ModDNA = MOD_DNA_STRING_CODEC,
         ModRNA = MOD_RNA_STRING_CODEC,
         NULL
  )
}

modsshortnames <- function(x)
{
  switch(x,
         ModDNA = .ModDNAorRNAcodes(MOD_DNA_BASE_CODES, c(), TRUE, TRUE, FALSE,
                                    col_name = "short_name"),
         ModRNA = .ModDNAorRNAcodes(MOD_RNA_BASE_CODES, c(), TRUE, TRUE, FALSE,
                                    col_name = "short_name"),
         NULL
  )
}

modsnomenclature <- function(x)
{
  switch(x,
         ModDNA = .ModDNAorRNAcodes(MOD_DNA_BASE_CODES, c(), TRUE, TRUE, FALSE,
                                    col_name = "nc"),
         ModRNA = .ModDNAorRNAcodes(MOD_RNA_BASE_CODES, c(), TRUE, TRUE, FALSE,
                                    col_name = "nc"),
         NULL
  )
}

comparable_seqtypes <- function(seqtype1, seqtype2)
{
  is_nucleo1 <- seqtype1 %in% c("ModDNA","ModRNA","DNA","RNA")
  is_nucleo2 <- seqtype2 %in% c("ModDNA","ModRNA","DNA","RNA")
  is_same_type <- (seqtype1 %in% c("ModDNA","DNA") & 
                     seqtype2 %in% c("ModDNA","DNA")) |
    (seqtype1 %in% c("ModRNA","RNA") & 
       seqtype2 %in% c("ModRNA","RNA"))
  is_nucleo1 == is_nucleo2 & is_same_type
}
