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
.normarg_nrec <- Biostrings:::.normarg_nrec
.normarg_skip <- Biostrings:::.normarg_skip
.normargCollapse <- Biostrings:::.normargCollapse
.normargLetters <- Biostrings:::.normargLetters
.normargOR <- Biostrings:::.normargOR
.normargWidth <- Biostrings:::.normargWidth
.normargWithIndels <- Biostrings:::normargWithIndels
.normargFixed <- Biostrings:::normargFixed
.normargAlgorithm <- Biostrings:::normargAlgorithm
.normargMaxMismatch <- Biostrings:::normargMaxMismatch
.normargMinMismatch <- Biostrings:::normargMinMismatch
.normargUseNames <- Biostrings:::normargUseNames
.normarg_at1 <- Biostrings:::.normarg_at1
.normarg_at2 <- Biostrings:::.normarg_at2

.close_filexp_list <- Biostrings:::.close_filexp_list
.to.ans.type <- Biostrings:::.to.ans.type
.xsbaseclass <- Biostrings:::xsbaseclass
.fromXStringViewsToStringSet <- Biostrings:::fromXStringViewsToStringSet
.isCharacterAlgo <- Biostrings:::isCharacterAlgo
.selectAlgo <- Biostrings:::selectAlgo
.character.matchPattern <- Biostrings:::.character.matchPattern
.compute_sorted_fasta_blocks_from_ssorted_fasta_index <- 
  Biostrings:::.compute_sorted_fasta_blocks_from_ssorted_fasta_index
.check_fasta_index <- Biostrings:::.check_fasta_index
.toSeqSnippet <- Biostrings:::toSeqSnippet
.XStringSet.show_frame_header <- Biostrings:::.XStringSet.show_frame_header
.V_recycle <- Biostrings:::.V_recycle
.H_recycle <- Biostrings:::.H_recycle
.make_IRanges_from_at <- Biostrings:::.make_IRanges_from_at

.DNA_BASE_CODES <- Biostrings:::DNA_BASE_CODES
.RNA_BASE_CODES <- Biostrings:::RNA_BASE_CODES
.XStringSet <- Biostrings:::XStringSet
.fromXStringViewsToStringSet <- Biostrings:::fromXStringViewsToStringSet
.XStringSetAsViews <- Biostrings:::.XStringSetAsViews
.XString.equal <- Biostrings:::.XString.equal
.XStringViews.equal <- Biostrings:::XStringViews.equal
.QualityScaledXStringSet <- Biostrings:::QualityScaledXStringSet

# IRanges
.new_Views <- IRanges:::new_Views
