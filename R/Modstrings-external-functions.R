#' @include Modstrings.R
NULL

# all these functions are used within Modstrings. However a reimplementation
# would at this time only cause code duplication. In addition this provides
# the benefit, that behaviour from the Biostrings package remains untouched.

XString <- Biostrings:::XString
XStringSet <- Biostrings:::XStringSet
XStringSetList <- Biostrings:::XStringSetList
.DNA_BASE_CODES <- Biostrings:::DNA_BASE_CODES
.RNA_BASE_CODES <- Biostrings:::RNA_BASE_CODES
.XStringSetAsViews <- Biostrings:::.XStringSetAsViews
.XString.equal <- Biostrings:::.XString.equal
.XStringViews.equal <- Biostrings:::XStringViews.equal
.QualityScaledXStringSet <- Biostrings:::QualityScaledXStringSet
.XString.nucleotide_frequency <- Biostrings:::.XString.nucleotide_frequency
.XStringSet.nucleotide_frequency <- 
  Biostrings:::.XStringSet.nucleotide_frequency

# import/export
.open_input_files <- XVector:::open_input_files
.normarg_compress <- XVector:::.normarg_compress
.normarg_nrec <- Biostrings:::.normarg_nrec
.normarg_skip <- Biostrings:::.normarg_skip
.close_filexp_list <- Biostrings:::.close_filexp_list
.compute_sorted_fasta_blocks_from_ssorted_fasta_index <- 
  Biostrings:::.compute_sorted_fasta_blocks_from_ssorted_fasta_index
.check_fasta_index <- Biostrings:::.check_fasta_index

# show functions
.namesW <- Biostrings:::.namesW
.toSeqSnippet <- Biostrings:::toSeqSnippet
.compact_ellipsis <- intToUtf8(0x2026L)
.XStringSet.show_frame_header <- Biostrings:::.XStringSet.show_frame_header

# IRanges
.new_Views <- IRanges:::new_Views
