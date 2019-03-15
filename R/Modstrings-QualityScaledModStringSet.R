#' @include Modstrings.R
#' @include Modstrings-ModStringSet.R
NULL

#' @name QualityScaledModStringSet
#' @aliases show,QualityScaledModStringSet-method
#' 
#' @title QualityScaledModDNAStringSet and QualityScaledModRNAStringSet objects
#' 
#' @description 
#' title
#' 
#' @param x For the \code{QualityScaled*StringSet} constructors: Either a 
#' character vector, or an \code{ModString}, \code{ModStringSet} or 
#' \code{ModStringViews} object.
#' 
#' For \code{writeQualityScaledXStringSet}: A 
#' \code{\link{QualityScaledModDNAStringSet}} or 
#' \code{\link{QualityScaledModRNAStringSet}} object.
#' @param quality A 
#' \code{\link[Biostrings:XStringQuality-class]{XStringQuality}} object.
#' @param filepath,nrec,skip,seek.first.rec,use.names,append,compress,compression_level
#'  See \code{\link[Biostrings:QualityScaledXStringSet-class]{QualityScaledXStringSet-class}}.
#' @param quality.scoring Specify the quality scoring used in the FASTQ file. 
#' Must be one of "phred" (the default), "solexa", or "illumina". If set to "
#' phred" (or "solexa" or "illumina"), the qualities will be stored in a 
#' \code{PhredQuality} (or \code{SolexaQuality} or \code{IlluminaQuality}, 
#' respectively) object.
#'  
#' @return a \code{QualityScaledModDNAStringSet} or 
#' \code{QualityScaledModDNAStringSet} object
NULL

# QualityScaledXStringSet objects

setClass("QualityScaledModStringSet",
         contains="QualityScaledXStringSet")

#' @rdname QualityScaledModStringSet
#' @export
setClass("QualityScaledModDNAStringSet",
         contains=c("ModDNAStringSet", "QualityScaledModStringSet")
)
#' @rdname QualityScaledModStringSet
#' @export
setClass("QualityScaledModRNAStringSet",
         contains=c("ModRNAStringSet", "QualityScaledModStringSet")
)
#' @rdname QualityScaledModStringSet
#' @export
QualityScaledModDNAStringSet <- function(x, quality) {
  .QualityScaledXStringSet(ModDNAStringSet(x), quality)
}
#' @rdname QualityScaledModStringSet
#' @export
QualityScaledModRNAStringSet <- function(x, quality){
  .QualityScaledXStringSet(ModRNAStringSet(x), quality)
}

setMethod("show", "QualityScaledModStringSet",
          function(object){
            cat("  A ", class(object), " instance containing:\n", sep="")
            cat("\n")
            selectMethod("show", "ModStringSet")(
              as(object, 
                 paste0(seqtype(object),"StringSet")))
            cat("\n")
            show(quality(object))
            cat("\n")
          }
)

# readQualityScaledDNAStringSet() / writeQualityScaledXStringSet() -------------

#' @rdname QualityScaledModStringSet
#' @export
readQualityScaledModDNAStringSet <- function(filepath,
                                             quality.scoring = c("phred",
                                                                 "solexa",
                                                                 "illumina"),
                                             nrec = -1L, skip = 0L,
                                             seek.first.rec = FALSE,
                                             use.names = TRUE)
{
  quality.scoring <- match.arg(quality.scoring)
  x <- readModDNAStringSet(filepath,
                           format = "fastq",
                           nrec,
                           skip,
                           seek.first.rec,
                           use.names,
                           with.qualities = TRUE)
  qualities <- mcols(x)[ , "qualities"]
  quals <- switch(quality.scoring,
                  phred = PhredQuality(qualities),
                  solexa = SolexaQuality(qualities),
                  illumina = IlluminaQuality(qualities))
  QualityScaledModDNAStringSet(x, quals)
}

#' @rdname QualityScaledModStringSet
#' @export
readQualityScaledModRNAStringSet <- function(filepath,
                                             quality.scoring = c("phred",
                                                                 "solexa",
                                                                 "illumina"),
                                             nrec = -1L,
                                             skip = 0L,
                                             seek.first.rec = FALSE,
                                             use.names = TRUE)
{
  quality.scoring <- match.arg(quality.scoring)
  x <- readModRNAStringSet(filepath,
                           format = "fastq",
                           nrec,
                           skip,
                           seek.first.rec,
                           use.names,
                           with.qualities = TRUE)
  qualities <- mcols(x)$qualities
  quals <- switch(quality.scoring,
                  phred = PhredQuality(qualities),
                  solexa = SolexaQuality(qualities),
                  illumina = IlluminaQuality(qualities))
  QualityScaledModRNAStringSet(x, quals)
}

#' @rdname QualityScaledModStringSet
#' @export
writeQualityScaledModStringSet <- function(x, filepath, append = FALSE,
                                           compress = FALSE,
                                           compression_level = NA)
{
  if (!is(x, "QualityScaledXStringSet")){
    stop("'x' must be a QualityScaledXStringSet object",
         call. = FALSE)
  }
  writeModStringSet(x,
                    filepath,
                    append,
                    compress,
                    compression_level,
                    format = "fastq",
                    qualities = quality(x))
}

