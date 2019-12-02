#' @include Modstrings.R
#' @include Modstrings-ModStringCodec.R
NULL

# These functions need to be here to access the modified functions of
# - .charToModString

#' @name ModString
#' @aliases ModString,MaskedModString-method ModString,AsIs-method 
#' ModString,ModString-method ModString,XString-method 
#' ModString,character-method ModString,factor-method
#' as.character,ModString-method as.vector,ModString-method
#' ==,ModString,ModString-method ==,ModString,XString-method
#' ==,XString,ModString-method
#' 
#' @title ModString objects
#' 
#' @description 
#' The virtual \code{ModString} class derives from the \code{XString} virtual
#' class. Like its parent and its children, it is used for storing sequences of
#' characters. However, the \code{XString}/\code{BString} class requires single
#' byte characters as the letters of the input sequences. The \code{ModString}
#' extends the capability for multi-byte chracters by encoding these characters
#' into a single byte characters using a dictionary for internal conversion. It
#' also takes care of different encoding behavior of operating systems.
#'
#' The \code{\link{ModDNAString}} and \code{\link{ModRNAString}} classes derive
#' from the \code{ModString} class and use the functionality to store nucleotide
#' sequences containing modified nucleotides. To describe modified RNA and DNA
#' nucleotides with a single letter, special characters are commonly used, eg.
#' from the greek alphabet, which are multi-byte characters.
#'
#' The \code{ModString} class is virtual and it cannot be directly used to
#' create an object. Please have a look at \code{\link{ModDNAString}} and
#' \code{\link{ModRNAString}} for the specific alphabets of the individual
#' classes.
NULL

#' @name ModDNAString
#' 
#' @title ModDNAString class
#' 
#' @description 
#' A \code{ModDNAString} object allows DNA sequences with modified nucleotides
#' to be stored and manipulated.
#' 
#' @details 
#' The ModDNAString class contains the virtual \code{\link{ModString}} class,
#' which is itself based on the \code{\link[Biostrings:XString-class]{XString}}
#' class. Therefore, functions for working with \code{XString} classes are 
#' inherited.
#' 
#' The \code{\link{alphabet}} of the ModDNAString class consist of the 
#' non-extended IUPAC codes "A,G,C,T,N", the gap letter "-", the hard masking 
#' letter "+", the not available letter "." and letters for individual 
#' modifications: \code{alphabet(ModDNAString())}.
#' 
#' Since the special characters are encoded differently depending on the OS and
#' encoding settings of the R session, it is not always possible to enter a DNA
#' sequence containing modified nucleotides via the R console. The most 
#' convinient solution for this problem is to use the function 
#' \code{\link{modifyNucleotides}} and modify and existing DNAString or
#' ModDNAString object.
#' 
#' A \code{ModDNAString} object can be converted into a \code{DNAString} object
#' using the \code{DNAstring()} constructor. Modified nucleotides are 
#' automaitcally converted intro their base nucleotides.
#' 
#' If a modified DNA nucleotide you want to work with is not part of the
#' alphabet, please let us know.
#'
#' @param x the input as a \code{character}.
#' @param start the postion in the character vector to use as start position in
#' the \code{ModDNAString} object (default \code{start = 1}).
#' @param nchar the width of the character vector to use in the 
#' \code{ModDNAString} object (default \code{nchar = NA}). The end position is
#' calculated as \code{start + nchar - 1}.
#' 
#' @return a \code{ModDNAString} object
#'
#' @examples
#' # Constructing ModDNAString containing an m6A
#' md1 <- ModDNAString("AGCT`")
#' 
#' # the alphabet of the ModDNAString class
#' alphabet(md1)
#' # due to encoding issues the shortNames can also be used
#' shortName(md1)
#' # due to encoding issues the nomenclature can also be used
#' nomenclature(md1) 
#' 
#' # convert to DNAString
#' d1 <- DNAString(md1)
NULL

#' @name ModRNAString
#' 
#' @title ModDNAString class
#' 
#' @description 
#' A \code{ModRNAString} object allows RNA sequences with modified nucleotides
#' to be stored and manipulated.
#' 
#' @details 
#' The ModRNAString class contains the virtual \code{\link{ModString}} class,
#' which is itself based on the \code{\link[Biostrings:XString-class]{XString}}
#' class. Therefore, functions for working with \code{XString} classes are 
#' inherited.
#' 
#' The alphabet of the ModRNAString class consist of the non-extended IUPAC 
#' codes "A,G,C,U", the gap letter "-", the hard masking letter "+", the not 
#' available letter "." and letters for individual modifications: 
#' \code{alphabet(ModRNAString())}.
#' 
#' Since the special characters are encoded differently depending on the OS and
#' encoding settings of the R session, it is not always possible to enter a RNA
#' sequence containing modified nucleotides via the R console. The most 
#' convinient solution for this problem is to use the function 
#' \code{\link{modifyNucleotides}} and modify and existing RNAString or
#' ModRNAString object.
#' 
#' A \code{ModRNAString} object can be converted into a \code{RNAString} object
#' using the \code{RNAstring()} constructor. Modified nucleotides are 
#' automaitcally converted intro their base nucleotides.
#' 
#' If a modified RNA nucleotide you want to work with is not part of the
#' alphabet, please let us know.
#'
#' @param x the input as a \code{character}.
#' @param start the postion in the character vector to use as start position in
#' the \code{ModRNAString} object (default \code{start = 1}).
#' @param nchar the width of the character vector to use in the 
#' \code{ModRNAString} object (default \code{nchar = NA}). The end position is
#' calculated as \code{start + nchar - 1}.
#' 
#' @return a \code{ModRNAString} object
#'
#' @examples
#' # Constructing ModDNAString containing an m6A and a dihydrouridine
#' mr1 <- ModRNAString("AGCU`D")
#' 
#' # the alphabet of the ModRNAString class
#' alphabet(mr1)
#' # due to encoding issues the shortNames can also be used
#' shortName(mr1)
#' # due to encoding issues the nomenclature can also be used
#' nomenclature(mr1)
#' 
#' # convert to RNAString
#' r1 <- RNAString(mr1)
NULL

# derived from Biostrings/R/XString-class.R ------------------------------------

setClass("ModString", contains = "XString")

#' @rdname ModString
#' @export
setClass("ModDNAString", contains = "ModString")
#' @rdname ModString
#' @export
setClass("ModRNAString", contains = "ModString")

#' @export
setMethod("seqtype", "ModDNAString", function(x) "ModDNA")
#' @export
setMethod("seqtype", "ModRNAString", function(x) "ModRNA")
#' @export
setReplaceMethod(
  "seqtype", "ModString",
  function(x, value)
  {
    ans_class <- paste0(value, "String")
    if(is(x,ans_class)){
      return(x)
    }
    ans_seq <- 
      XVector:::extract_character_from_XRaw_by_ranges(x, 1L, length(x),
                                                      collapse=FALSE,
                                                      lkup=NULL)
    ans_seq <-
      .convert_one_byte_codes_to_originating_base(ans_seq,
                                                  modscodec(seqtype(x)))
    do.call(ans_class,list(ans_seq))
  }
)

# low level functions construct XString objects and extract character

setMethod(
  "extract_character_from_XString_by_positions", "ModString",
  function(x, pos, collapse=FALSE)
  {
    ans <- callNextMethod()
    codec <- modscodec(seqtype(x))
    .convert_one_byte_codes_to_letters(ans, codec)
  }
)

setMethod(
  "extract_character_from_XString_by_ranges", "ModString",
  function(x, start, width, collapse=FALSE)
  {
    ans <- callNextMethod()
    codec <- modscodec(seqtype(x))
    .convert_one_byte_codes_to_letters(ans, codec)
  }
)

setMethod(
  "make_XString_from_string", "ModString",
  function(x0, string, start, width)
  {
    codec <- modscodec(seqtype(x0))
    string <- .convert_letters_to_one_byte_codes(string, codec)
    callNextMethod()
  }
)

# Constructor ------------------------------------------------------------------

#' @rdname ModDNAString
#' @export
ModDNAString <- function(x = "", start = 1, nchar = NA){
  XString("ModDNA", x, start = start, width = nchar)
}
#' @rdname ModRNAString
#' @export
ModRNAString <- function(x = "", start = 1, nchar = NA){
  XString("ModRNA", x, start = start, width = nchar)
}

# Coercion ---------------------------------------------------------------------

#' @rdname Modstrings-internals
#' @export
setMethod(
  "XString", "ModString",
  function(seqtype, x, start = NA, end = NA, width = NA)
  {
    ans <- subseq(x, start = start, end = end, width = width)
    ans_class <- paste0(seqtype, "String")
    if(is(ans,ans_class)){
      return(ans)
    }
    # convert over "base" classes to convert T/U
    seqtype(ans) <- gsub("Mod","",seqtype(ans))
    seqtype(ans) <- gsub("Mod","",seqtype)
    seqtype(ans) <- seqtype
    ans
  }
)

#' @export
setAs("XString", "ModDNAString",
      function(from) {
        seqtype(from) <- "ModDNA"
        from
      }
)
#' @export
setAs("XString", "ModRNAString",
      function(from) {
        seqtype(from) <- "ModRNA"
        from
      }
)

#' @export
setAs("character", "ModDNAString", function(from) ModDNAString(from))
#' @export
setAs("character", "ModRNAString", function(from) ModRNAString(from))
#' @export
setMethod(
  "as.vector", "ModString",
  function(x){
    codec <- modscodec(seqtype(x))
    x_alphabet <- letters(codec)
    code2pos <- as.integer(unlist(lapply(oneByteCodes(codec), charToRaw)))
    x_alphabet <- x_alphabet[order(code2pos)]
    ans <- as.integer(x)
    attributes(ans) <- list(levels = x_alphabet, class = "factor")
    as.vector(ans)
  }
)

# Show 

### Placeholder, initialized in .onLoad()
MODDNA_COLORED_LETTERS <- NULL
MODRNA_COLORED_LETTERS <- NULL

### Return a named character vector where all the names are single letters.
### Colors for A, C, G, and T were inspired by
###   https://en.wikipedia.org/wiki/Nucleotide#Structure
### Called in .onLoad() to initialize MODDNA_COLORED_LETTERS and 
### MODRNA_COLORED_LETTERS
#' @importFrom crayon make_style inverse
make_MODDNA_COLORED_LETTERS <- function()
{
  ## modified DNA nucleotide letters
  ans <- setNames(sprintf(make_style(rgb(0.2,0.2,0.2), bg=TRUE)(make_style(rgb(0,1,0))("%s")),
                          letters(MOD_DNA_STRING_CODEC)),
                  letters(MOD_DNA_STRING_CODEC))
  # base colours
  ans["A"] <- make_style(rgb(1, 0.5, 0.5), bg=TRUE)(make_style("black")("A"))
  ans["C"] <- make_style(rgb(0.5, 1, 0.5), bg=TRUE)(make_style("black")("C"))
  ans["G"] <- make_style(rgb(0.5, 1, 1), bg=TRUE)(make_style("black")("G"))
  ans["T"] <- make_style(rgb(1, 0.8, 0.5), bg=TRUE)(make_style("black")("T"))
  ans["N"] <- make_style("grey", bg=TRUE)(make_style(rgb(0.4,0.4,0.4))("N"))
  ans["-"] <- "-"
  ans["+"] <- "+"
  ans["."] <- "."
  ans
}
make_MODRNA_COLORED_LETTERS <- function()
{
  ## modified RNA nucleotide letters
  ans <- setNames(sprintf(make_style(rgb(0.2,0.2,0.2), bg=TRUE)(make_style(rgb(0,1,0))("%s")),
                          letters(MOD_RNA_STRING_CODEC)),
                  letters(MOD_RNA_STRING_CODEC))
  # base colours
  ans["A"] <- make_style(rgb(1, 0.5, 0.5), bg=TRUE)(make_style("black")("A"))
  ans["C"] <- make_style(rgb(0.5, 1, 0.5), bg=TRUE)(make_style("black")("C"))
  ans["G"] <- make_style(rgb(0.5, 1, 1), bg=TRUE)(make_style("black")("G"))
  ans["U"] <- make_style(rgb(1, 0.8, 0.5), bg=TRUE)(make_style("black")("U"))
  ans["N"] <- make_style("grey", bg=TRUE)(make_style(rgb(0.4,0.4,0.4))("N"))
  ans["-"] <- "-"
  ans["+"] <- "+"
  ans["."] <- "."
  ans
}

### 'x' must be a character vector.
.add_modx_colors <- function(x,COLORED_LETTERS){
  ans <- vapply(x,
                function(xi){
                  xi <- strsplit(xi,"")[[1L]]
                  m <- match(xi, names(COLORED_LETTERS))
                  match_idx <- which(!is.na(m))
                  xi[match_idx] <- COLORED_LETTERS[m[match_idx]]
                  paste0(xi, collapse="")
                },
                character(1),
                USE.NAMES=FALSE
  )
  x_names <- names(x)
  if (!is.null(x_names))
    names(ans) <- x_names
  ans
}
.add_moddna_colors <- function(x)
{
  .add_modx_colors(x, MODDNA_COLORED_LETTERS)
}

.add_modrna_colors <- function(x)
{
  .add_modx_colors(x, MODRNA_COLORED_LETTERS)
}

add_colors <- function(x) UseMethod("add_colors")
add_colors.default <- identity
add_colors.ModDNA <- .add_moddna_colors
add_colors.ModRNA <- .add_modrna_colors

setMethod("show", "ModString",
          function(object)
          {
            object_len <- object@length
            cat(object_len, "-letter ", class(object), " object\n", sep="")
            snippet <- .toSeqSnippet(object, getOption("width") - 5L)
            cat("seq: ", add_colors(snippet), "\n", sep="")
          }
)

# Comparison -------------------------------------------------------------------

.compare_ModString <- function(e1,
                               e2){
  if (!comparable_seqtypes(seqtype(e1), seqtype(e2))) {
    class1 <- class(e1)
    class2 <- class(e2)
    stop("comparison between a \"", class1, "\" instance ",
         "and a \"", class2, "\" instance ",
         "is not supported")
  }
  if(!is(e1,"ModString")){
    e1 <- BString(e1)
  }
  if(!is(e2,"ModString")){
    e2 <- BString(e2)
  }
  .XString.equal(e1, e2)
}

#' @export
setMethod("==", signature(e1 = "ModString", e2 = "ModString"),
          function(e1, e2) .compare_ModString(e1, e2)
)
#' @export
setMethod("==", signature(e1 = "ModString", e2 = "XString"),
          function(e1, e2) .compare_ModString(e1, e2)
)
#' @export
setMethod("==", signature(e1 = "XString", e2 = "ModString"),
          function(e1, e2) .compare_ModString(e1, e2)
)

# these accessors are not provided by the XVector package
setGeneric(name = "sharedXVector",
           signature = "x",
           def = function(x) standardGeneric("sharedXVector"))
setGeneric(name = "offsetXVector",
           signature = "x",
           def = function(x) standardGeneric("offsetXVector"))

setMethod("sharedXVector","ModString",
          function(x) x@shared)
setMethod("offsetXVector","ModString",
          function(x) x@offset)
