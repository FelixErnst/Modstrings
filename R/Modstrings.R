#' @title Modstrings: implementation of Biostrings to work with nucleotide 
#' sequences containing modified nucleotides.
#'
#' @author Felix G M Ernst [aut,cre] and Denis L.J. Lafontaine [ctb]
#'
#' @description
#' Representing nucleotide modifications in a nucleotide sequence is usually
#' done via special characters from a number of sources. This represents a
#' challenge to work with in R and the \code{Biostrings} package. The
#' \code{Modstrings} package implements this functionallity for RNA and DNA
#' sequences containing modified nucleotides by translating the character
#' internally in order to work with the infrastructure of the \code{Biostrings}
#' package. For this the \code{ModRNAString} and \code{ModDNAString} classes and
#' derivates and functions to construct and modify these objects despite the
#' encoding issues are implemenented. In addition the conversion from sequences
#' to list like location information (and the reverse operation) is implemented
#' as well.
#' 
#' A good place to start would be the vignette and the man page for the
#' \code{\link{ModStringSet}} objects.
#' 
#' The alphabets for the modifications used in this package are based on the 
#' compilation of RNA modifications by \url{http://modomics.genesilico.pl} by
#' the Bujnicki lab and DNA modifications \url{https://dnamod.hoffmanlab.org}
#' by the Hoffman lab. Both alphabets were modified to remove some incompatible
#' characters.
#'
#' @docType package
#' @name Modstrings
NULL

#' @import methods
#' @import BiocGenerics
#' @import Biostrings
#' @import GenomicRanges
#' @import S4Vectors
NULL

#' @name Modstrings-internals
#' @aliases seqtype,ModDNAString-method seqtype,ModRNAString-method 
#' seqtype<-,ModString-method seqtype<-,ModStringSet-method
#' 
#' @title Modstrings internals
#' 
#' @param seqtype,x,start,end,width,use.names used internally
#' 
#' @return a XString* object
#' 
#' @keywords internal
#' 
#' @description 
#' Analog to \code{Biostrings} there are a few functions, which should only
#' be used internally. Otherwise take care.
NULL
