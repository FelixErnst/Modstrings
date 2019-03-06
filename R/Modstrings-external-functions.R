#' @include Modstrings.R
NULL

# all these functions are used within Modstrings. However a reimplementation
# would at this time only cause code duplication. In addition this provides
# the benefit, that behaviour from the Biostrings package remains untouched.

# XVectors 
.open_input_files <- XVector:::open_input_files
.normarg_compress <- XVector:::.normarg_compress
.new_XVectorList_from_list_of_XVector <-
  XVector:::new_XVectorList_from_list_of_XVector

# Biostrings
.DNA_BASE_CODES <- Biostrings:::DNA_BASE_CODES
.RNA_BASE_CODES <- Biostrings:::RNA_BASE_CODES
.XString <- Biostrings:::XString
.XStringSet <- Biostrings:::XStringSet
.XStringSetAsViews <- Biostrings:::.XStringSetAsViews
.XString.equal <- Biostrings:::.XString.equal
.XStringViews.equal <- Biostrings:::XStringViews.equal
.QualityScaledXStringSet <- Biostrings:::QualityScaledXStringSet
.fromXStringViewsToStringSet <- Biostrings:::fromXStringViewsToStringSet

# import/export
.normarg_nrec <- Biostrings:::.normarg_nrec
.normarg_skip <- Biostrings:::.normarg_skip
.close_filexp_list <- Biostrings:::.close_filexp_list
.compute_sorted_fasta_blocks_from_ssorted_fasta_index <- 
  Biostrings:::.compute_sorted_fasta_blocks_from_ssorted_fasta_index
.check_fasta_index <- Biostrings:::.check_fasta_index

# show functions
.toSeqSnippet <- Biostrings:::toSeqSnippet
.XStringSet.show_frame_header <- Biostrings:::.XStringSet.show_frame_header

# tied to C calls imported from Biostrings
.normargCollapse <- Biostrings:::.normargCollapse
.normargLetters <- Biostrings:::.normargLetters
.normargOR <- Biostrings:::.normargOR
.normargWidth <- Biostrings:::.normargWidth
.normargWithIndels <- Biostrings:::normargWithIndels
.normargSubject <- Biostrings:::normargSubject
.normargFixed <- Biostrings:::normargFixed
.normargAlgorithm <- Biostrings:::normargAlgorithm
.normargMaxMismatch <- Biostrings:::normargMaxMismatch
.normargMinMismatch <- Biostrings:::normargMinMismatch
.normargUseNames <- Biostrings:::normargUseNames
.to.ans.type <- Biostrings:::.to.ans.type
.xsbaseclass <- Biostrings:::xsbaseclass
.isCharacterAlgo <- Biostrings:::isCharacterAlgo
.selectAlgo <- Biostrings:::selectAlgo
.character.matchPattern <- Biostrings:::.character.matchPattern

# IRanges
.new_Views <- IRanges:::new_Views
