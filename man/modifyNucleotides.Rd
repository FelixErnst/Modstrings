% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/Modstrings-modifyNucleotide.R
\name{modifyNucleotides}
\alias{modifyNucleotides}
\alias{modifyNucleotides,ModString-method}
\alias{modifyNucleotides,ModStringSet-method}
\alias{modifyNucleotides,DNAString-method}
\alias{modifyNucleotides,RNAString-method}
\alias{modifyNucleotides,DNAStringSet-method}
\alias{modifyNucleotides,RNAStringSet-method}
\title{Modifying nucleotides in a nucleotide sequence (or set of sequences) 
at specified locations}
\usage{
modifyNucleotides(
  x,
  at,
  mod,
  nc.type = "short",
  stop.on.error = TRUE,
  verbose = FALSE
)

\S4method{modifyNucleotides}{ModString}(
  x,
  at,
  mod,
  nc.type = c("short", "nc"),
  stop.on.error = TRUE,
  verbose = FALSE
)

\S4method{modifyNucleotides}{ModStringSet}(
  x,
  at,
  mod,
  nc.type = c("short", "nc"),
  stop.on.error = TRUE,
  verbose = FALSE
)

\S4method{modifyNucleotides}{DNAString}(
  x,
  at,
  mod,
  nc.type = c("short", "nc"),
  stop.on.error = TRUE,
  verbose = FALSE
)

\S4method{modifyNucleotides}{RNAString}(
  x,
  at,
  mod,
  nc.type = c("short", "nc"),
  stop.on.error = TRUE,
  verbose = FALSE
)

\S4method{modifyNucleotides}{DNAStringSet}(
  x,
  at,
  mod,
  nc.type = c("short", "nc"),
  stop.on.error = TRUE,
  verbose = FALSE
)

\S4method{modifyNucleotides}{RNAStringSet}(
  x,
  at,
  mod,
  nc.type = c("short", "nc"),
  stop.on.error = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{x}{a \code{\link{ModString}} or \code{\link{ModStringSet}} object}

\item{at}{the location where the modification should be made.

The same input as in the original \code{\link[Biostrings]{replaceLetterAt}}
are expected:

If x is a \code{\link{ModString}} object, then at is typically an integer
vector with no NAs but a logical vector or Rle object is valid too. Locations
can be repeated and in this case the last replacement to occur at a given
location prevails.

If x is a rectangular \code{\link{ModStringSet}} object, then \code{at} must
be a matrix of logicals with the same dimensions as x. If the
\code{\link{ModStringSet}} is not rectangular, \code{at} must be a list of
logical vectors.}

\item{mod}{The modification short name or nomenclature

If \code{x} is a \code{\link{ModString}} object, then letter must be a 
\code{\link{ModString}} object or a character vector (with no \code{NA}) with
a total number of letters \code{(sum(nchar(letter)))} equal to the number of
locations specified in at.

If \code{x} is a rectangular \code{\link{ModStringSet}} object, then letter
must be a \code{\link{ModStringSet}} object, a list of character vectors or a
\code{\link[IRanges:AtomicList-class]{CharacterList}} of the same length as
x. In addition, the number of letters in each element of letter must match
the number of locations specified in the corresponding row of at
\code{(all(width(letter) == rowSums(at)))}.}

\item{nc.type}{the type of nomenclature to be used. Either "short" or "nc".
"Short" for m3C would be "m3C", "nc" for m3C would be "3C". (
\code{default = "short"})}

\item{stop.on.error}{For \code{combineIntoModstrings}: \code{TRUE}(default)
or \code{FALSE}: Should an error be raised upon encounter of incompatible 
positions?}

\item{verbose}{See \code{\link[Biostrings]{replaceLetterAt}}.}
}
\value{
the input \code{\link{ModString}} or \code{\link{ModStringSet}}
object with the changes applied
}
\description{
\code{\link{modifyNucleotides}} modifies a nucleotide in a sequence (or set 
of sequences) based on the type of modification provided. It checks for the 
identity of the base nucleotide to be
}
\examples{
# modify nucleotides in a ModDNAString 
seq <- ModDNAString("AGTC")
seq

mseq1 <- modifyNucleotides(seq,c(1,2,4),c("1mA","7mG","3mC"))
mseq1

# This fails since m7G requires a G at the selected position in the sequence
\dontrun{
mseq <- modifyNucleotides(seq,c(3),c("7mG"))
}

# modify nucleotides in a ModRNAString 
seq <- ModRNAString("AGUC")
seq

mseq1 <- modifyNucleotides(seq,c(1,2,4),c("m1A","m7G","m3C"))
mseq1

# This fails since m7G requires a G at the selected position in the sequence
\dontrun{
mseq <- modifyNucleotides(seq,c(3),c("m7G"))
}
}
