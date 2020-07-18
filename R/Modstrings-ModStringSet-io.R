#' @include Modstrings.R
#' @include Modstrings-ModString.R
NULL

#' @name ModStringSet-io
#' @aliases readModDNAStringSet readModRNAStringSet writeModStringSet
#' 
#' @title Read/write an ModStringSet object from/to a file
#' 
#' @description 
#' Functions to read/write an ModStringSet object from/to a file.
#'
#' @param x A \code{ModStringSet} object.
#' @param filepath,format,nrec,skip,seek.first.rec,use.names,with.qualities,append,compress,compression_level,...
#'  See \code{\link[Biostrings:XStringSet-io]{XStringSet-io}} for more details.
#' 
#' @return A \code{\link{ModStringSet}} of the defined type.
#'
#' @examples
#' seqs <- paste0(paste(alphabet(ModDNAString()), collapse = ""),
#'                c("A","G","T"))
#' set <- ModDNAStringSet(seqs)
#' file <- tempfile()
#' writeModStringSet(set, file)
#' read <- readModDNAStringSet(file)
NULL

# derived from Biostrings/R/XStringSet-io.R ------------------------------------

# fasta ------------------------------------------------------------------------

.fasta.index <- function(filepath,
                         nrec = -1L,
                         skip = 0L,
                         seek.first.rec = FALSE)
{
  filexp_list <- .open_input_files(filepath)
  on.exit(.close_filexp_list(filexp_list))
  nrec <- .normarg_nrec(nrec)
  skip <- .normarg_skip(skip)
  if (!.is_a_bool(seek.first.rec)){
    stop("'seek.first.rec' must be TRUE or FALSE", call. = FALSE)
  }
  # Guess what!? it does not like invalid one-letter sequence codes
  # suppress the warning. Results are valid noneheless for this purpose
  # However the offset cannot be used, since readBin output cannot be encoded
  # correctly
  ans <- .call_fasta_index(filexp_list, nrec, skip, seek.first.rec)
  ## 'expath' will usually be the same as 'filepath', except when 'filepath'
  ## contains URLs which will be replaced by the path to the downloaded file.
  expath <- vapply(filexp_list,
                   attr,
                   character(1),
                   "expath",
                   USE.NAMES=FALSE)
  ans$filepath <- expath[ans[ , "fileno"]]
  ans
}

# read a fasta sequence from a block of lines
.parse_fasta_block <- function(lines, seqtype)
{
  ans <- do.call(seqtype, list(paste(lines, collapse = "")))
  ans
}

#' @importFrom stringr str_trim str_sub
.read_ModStringSet_from_fasta_file_data <- function(con, seqtype,
                                                    nLines = 10000L)
{
  ans <- list()
  lines <- readLines(con,
                     nLines,
                     encoding = "UTF-8") 
  linesRead <- length(lines)
  if(linesRead== 0){
    warning("Empty file.", call. = FALSE)
  }
  # while is used, since the offset fromt he fai cannot be used to direct
  # readBin, since readBin result cannot be encoded correctly
  while(linesRead > 0){
    mergePrevious <- FALSE
    names <- NULL
    # get lines numbers for blocks of lines with one fasta sequence
    cuts <- unlist(lapply(seq_along(lines), function(i){
      if(stringr::str_sub(lines[i],1,1) == ">") return(i)
    }))
    # if no fasta block is detected, all lines read belong to the previous one
    if(!is.null(cuts)){
      # extract the names
      names <- stringr::str_trim(stringr::str_sub(lines[cuts],
                                                  2,
                                                  nchar(lines[cuts])))
      lines[cuts] <- ""
      # if a fasta block does not start at the beginning, the remaining seq 
      # belongs to the previous block
      if(cuts[1] == 1){
        cutsBegin <- cuts
      } else {
        cutsBegin <- c(1,cuts)
        mergePrevious <- TRUE
      }
    } else {
      cutsBegin <- 1
      mergePrevious <- TRUE
    }
    cutsEnd <- cutsBegin - 1
    cutsEnd <- c(cutsEnd[cutsEnd > 0],linesRead)
    # split lines into blocks and get fasta sequence
    f <- unlist(mapply(rep, 
                       cutsBegin, 
                       (cutsEnd - cutsBegin + 1), 
                       SIMPLIFY = FALSE))
    # get a list if seqtype strings
    ans_add <- 
      lapply(split(lines, f),
             function(nextLines){
               .parse_fasta_block(nextLines[stringr::str_trim(nextLines) != ""],
                                  seqtype)
             })
    # if the last of the already read in and first of the new ones need merging 
    if(mergePrevious){
      ans[[length(ans)]] <- Biostrings::xscat(ans[[length(ans)]],
                                              ans_add[[1]])
      ans_add[[1]] <- NULL
    }
    #
    names(ans_add) <- names
    ans <- append(ans,
                  ans_add)
    # load new lines if available
    lines <- readLines(con,
                       nLines,
                       encoding = "UTF-8")
    linesRead <- length(lines)
  }
  ans <- ans[!vapply(ans,is.null,logical(1))]
  ans
}

### Fasta index 'ssorted_fai' must be strictly sorted by "recno". This is NOT
### checked!
.read_ModStringSet_from_ssorted_fasta_index <- function(ssorted_fai,
                                                        elementType)
{
  ## Prepare 'nrec_list' and 'offset_list'.
  fasta_blocks <-
    .compute_sorted_fasta_blocks_from_ssorted_fasta_index(ssorted_fai)
  nrec_list <- split(fasta_blocks[ , "nrec"],
                     fasta_blocks[ , "fileno"],
                     drop = TRUE)
  ## Prepare 'used_filepath'.
  filepath <- ssorted_fai[ , "filepath"]
  fileno <- ssorted_fai[ , "fileno"]
  recno <- ssorted_fai[ , "recno"]
  used_fileno <- as.integer(names(nrec_list))
  used_filepath <- filepath[match(used_fileno, fileno)]
  # open files
  cons <- lapply(used_filepath, file, "rt")
  on.exit(lapply(cons,close))
  # parse files for ModStrings
  ans <- mapply(.read_ModStringSet_from_fasta_file_data,
                cons,
                MoreArgs = list(seqtype = elementType),
                SIMPLIFY = FALSE)
  # construct a ModStringSet
  ans <- ans[!vapply(ans,is.null,logical(1))]
  ans <- lapply(ans,
                function(seqs){
                  do.call(paste0(elementType,"Set"),list(seqs))
                })
  if(length(ans) == 0 ) {
    stop("No fasta blocks found. ",
         call. = FALSE)
  }
  ans <- do.call(c,ans)
  ans <- ans[recno]
  # check nrec against results and apply offset
  # this checks the found results against the number of .fasta.index
  nrecFound <- lapply(split(ans,fileno), length)
  if(any(unlist(nrecFound) != unlist(nrec_list))){
    warning("Not all fasta blocks were recognized. ",
            call. = FALSE)
  }
  ans
}

.read_ModStringSet_from_fasta_index <- function(fai, use.names, elementType)
{
  .check_fasta_index(fai)
  ## Create a "strictly sorted" version of 'fai' by removing duplicated rows
  ## and sorting the remaining rows by ascending "recno".
  recno <- fai[ , "recno"]
  ssorted_recno <- sort(unique(recno))
  ssorted_fai <- fai[match(ssorted_recno, recno), , drop = FALSE]
  ans <- .read_ModStringSet_from_ssorted_fasta_index(ssorted_fai, elementType)
  ## Re-order XStringSet object to make it parallel to 'recno'.
  ans <- ans[match(recno, ssorted_recno)]
  ## Set names on XStringSet object.
  if (use.names){
    names(ans) <- fai[ , "desc"]
  }
  ans
}

# Fastq ------------------------------------------------------------------------

.check_fastq_markup <- function(lines, nlines)
{
  if(nlines != 4){
    stop("Partial fastq block found. Aborting...")
  }
  if(substr(lines[1],1,1) != "@"){
    stop("Wrong format. '@' not found. Aborting...")
  }
  if(substr(lines[1],1,1) != "@"){
    stop("Wrong format. '@' not found. Aborting...")
  }
}

.get_names_from_fastq_block <- function(lines, with.qualities)
{
  name_seq <- substr(lines[1],2,nchar(lines[1]))
  if(with.qualities){
    name_qual <- substr(lines[3],2,nchar(lines[3]))
    if(name_seq != name_qual){
      warning("Name of sequence and quality do not match. '",
              name_seq,
              "'")
    }
  }
  name_seq
}

.get_sequences_from_fastq_block <- function(lines, seqtype, with.qualities)
{
  seq <- do.call(seqtype, list(lines[2]))
  if(with.qualities){
    attr(seq,"qualities") <- BString(lines[4])
  }
  seq
}

.read_ModStringSet_from_fastq_file_data <- function(con, seqtype,
                                                    with.qualities,
                                                    seek.first.rec, use.names)
{
  ans <- list()
  names <- NULL
  # Search for the first @
  if(seek.first.rec){
    lines <- ""
    while(substr(lines,1,1) != "@" & length(lines) != 0L){
      lines <- readLines(con,
                         1,
                         encoding = "UTF-8")
    }
    if(length(lines) == 0L){
      stop("No start markaup '@' for any fastq block found.", call. = FALSE)
    }
    lines <- c(lines,
               readLines(con,
                         3,
                         encoding = "UTF-8") )
  } else {
    lines <- readLines(con,
                       4,
                       encoding = "UTF-8")
  }
  linesRead <- length(lines)
  while(linesRead == 4){
    .check_fastq_markup(lines, linesRead)
    names <- c(names,.get_names_from_fastq_block(lines,
                                                 with.qualities))
    ans <- c(ans,
             .get_sequences_from_fastq_block(lines,
                                             seqtype,
                                             with.qualities))
    lines <- readLines(con,
                       4,
                       encoding = "UTF-8") 
    linesRead <- length(lines)
  }
  if(use.names){
    names(ans) <- names
  }
  ans <- ans[!vapply(ans,is.null,logical(1))]
  ans
}

.read_ModStringSet_from_fastq <- function(filepath, nrec, skip, seek.first.rec,
                                          use.names, elementType,
                                          with.qualities)
{
  cons <- lapply(filepath, file, "rt")
  on.exit(lapply(cons,close))
  nrec <- .normarg_nrec(nrec)
  skip <- .normarg_skip(skip)
  if (!.is_a_bool(seek.first.rec)){
    stop("'seek.first.rec' must be TRUE or FALSE", call. = FALSE)
  }
  if (!.is_a_bool(with.qualities)){
    stop("'with.qualities' must be TRUE or FALSE", call. = FALSE)
  }
  # parse files for ModStrings
  ans <- mapply(.read_ModStringSet_from_fastq_file_data,
                cons,
                MoreArgs = list(seqtype = elementType,
                                with.qualities = with.qualities,
                                seek.first.rec = seek.first.rec,
                                use.names = use.names),
                SIMPLIFY = FALSE)
  # construct a ModStringSet
  ans <- ans[!vapply(ans,is.null,logical(1))]
  if(length(ans) == 0 ) {
    stop("No fastq blocks found. ",
         call. = FALSE)
  }
  if(with.qualities){
    qual <- lapply(ans,
                   function(a){
                     unname(BStringSet(lapply(a, attr, "qualities")))
                   })
    ans <- do.call(c,lapply(ans,paste0(elementType,"Set")))
    qual <- do.call(c,qual)
    #
    f <- seq_along(ans)
    if(nrec > -1L){
      f <- f[f %in% (skip + 1L):(skip + 1L + nrec - 1L)]
    }
    ans <- ans[f]
    qual <- qual[f]
    mcols(ans)$qualities <- qual
    return(ans)
  }
  ans <- do.call(c,lapply(ans,paste0(elementType,"Set")))
  f <- seq_along(ans)
  if(nrec > -1L){
    f <- f[f %in% (skip + 1L):(skip + 1L + nrec - 1L)]
  }
  ans <- ans[f]
  ans
}

.read_ModStringSet <- function(filepath, format, nrec = -1L, skip = 0L,
                               seek.first.rec = FALSE, use.names = TRUE,
                               seqtype = "B", with.qualities = FALSE)
{
  if (!.is_a_string(format)){
    stop("'format' must be a single string", call. = FALSE)
  }
  format <- match.arg(tolower(format), c("fasta","fastq"))
  if (!.is_a_bool(use.names)){
    stop("'use.names' must be TRUE or FALSE", call. = FALSE)
  }
  elementType <- paste0(seqtype, "String")
  # fastq
  if (format == "fastq") {
    ans <- .read_ModStringSet_from_fastq(filepath, nrec, skip, seek.first.rec,
                                         use.names, elementType, with.qualities)
    return(ans)
  }
  # fasta
  if (!identical(with.qualities, FALSE)){
    stop("The 'with.qualities' argument is only supported ",
         "when reading a FASTQ file.",
         call. = FALSE)
  }
  if (is.data.frame(filepath)) {
    if (!(identical(nrec, -1L) &&
          identical(skip, 0L) &&
          identical(seek.first.rec, FALSE))) {
      warning("'nrec', 'skip', and 'seek.first.rec' are ",
              "ignored when 'filepath' is a data frame")
    }
    fai <- filepath
  } else {
    fai <- .fasta.index(filepath, nrec = nrec, skip = skip,
                        seek.first.rec = seek.first.rec)
  }
  .read_ModStringSet_from_fasta_index(fai, use.names, elementType)
}


#' @rdname ModStringSet-io
#' @export
readModDNAStringSet <- function(filepath, format = "fasta", nrec = -1L,
                                skip = 0L, seek.first.rec = FALSE,
                                use.names = TRUE, with.qualities = FALSE) {
  ans <- .read_ModStringSet(filepath, format = format, nrec = nrec, skip = skip,
                            seek.first.rec = seek.first.rec,
                            use.names = use.names, seqtype = "ModDNA",
                            with.qualities = with.qualities)
  ans
}

#' @rdname ModStringSet-io
#' @export
readModRNAStringSet <- function(filepath, format = "fasta", nrec = -1L,
                                skip = 0L, seek.first.rec = FALSE,
                                use.names = TRUE, with.qualities = FALSE) {
  ans <- .read_ModStringSet(filepath, format = format, nrec = nrec, skip = skip,
                            seek.first.rec = seek.first.rec,
                            use.names = use.names, seqtype = "ModRNA",
                            with.qualities = with.qualities)
  ans
}

# derived from Biostrings/R/XStringSet-io.R ------------------------------------
# Writing ModStringSets to file
# The R implementation is a proof of concept for now. C implementation would of
# course be faster. This would mean that the one byte lookup would be also done
# in C.

.format_as_fasta <- function(x, width){
  names <- paste0(">",names(x))
  x <- as.list(as.character(x))
  x <- mapply(function(z,
                       name){
                ans <- stringr::str_sub(z,
                                        start = seq.int(0,
                                                        nchar(z),
                                                        width) + 1,
                                        end = c(c(seq.int(0,
                                                          nchar(z),
                                                          width))[-1],
                                                nchar(z)))
                c(name,ans)
              },
              x,
              names,
              SIMPLIFY = FALSE)
  unlist(x)
}

.write_ModStringSet_to_fasta <- function(x, con, width = 80L, ...)
{
  if (!.is_a_number(width)){
    stop("'width' must be a single integer", call. = FALSE)
  }
  if (width < 1L){
    stop("'width' must be an integer >= 1", call. = FALSE)
  }
  length <- length(x)
  chunk_size <- 100
  no_of_chunks <- (length %/% chunk_size) + 1
  f <- rep(seq_len(no_of_chunks),
           c(rep(chunk_size,no_of_chunks - 1),length %% chunk_size))
  splitX <- split(x,f)
  lapply(splitX,
         function(z){
           str <- .format_as_fasta(z, width)
           writeLines(str,
                      con = con,
                      useBytes = TRUE)
         })
  NULL
}

.format_as_fastq <- function(x, qualities)
{
  namesSeq <- paste0("@",names(x))
  namesQual <- paste0("+",names(x))
  x <- as.list(as.character(x))
  qualities <- as.list(as.character(qualities))
  x <- mapply(function(z,
                       q,
                       nameZ,
                       nameQ){
    c(nameZ,z,nameQ,q)
  },
  x,
  qualities,
  namesSeq,
  namesQual,
  SIMPLIFY = FALSE)
  unlist(x)
}

.write_ModStringSet_to_fastq <- function(x, con, qualities = NULL, ...)
{
  if (is.null(qualities))
    qualities <- mcols(x)$qualities
  if (!is.null(qualities)) {
    if (!is(qualities, "BStringSet"))
      stop("'qualities' must be NULL or a BStringSet object",
           call. = FALSE)
    if (length(qualities) != length(x))
      stop("'x' and 'qualities' must have the same length",
           call. = FALSE)
  }
  length <- length(x)
  chunk_size <- 100
  no_of_chunks <- (length %/% chunk_size) + 1
  f <- rep(seq_len(no_of_chunks),
           c(rep(chunk_size,no_of_chunks - 1),length %% chunk_size))
  splitX <- split(x,f)
  if (is.null(qualities)) {
    qualities <- unname(BStringSet(unlist(lapply(mapply(rep,
                                                        "!",
                                                        width(x),
                                                        SIMPLIFY = FALSE),
                                                 paste,
                                                 collapse = ""))))
  }
  splitQualities <- split(qualities,f)
  mapply(
    function(z,q){
      str <- .format_as_fastq(z, q)
      writeLines(str,
                 con = con,
                 useBytes = TRUE)
    },
    splitX,
    splitQualities)
  NULL
}

#' @rdname ModStringSet-io
#' @export
writeModStringSet <- function(x, filepath, append = FALSE, compress = FALSE,
                              compression_level = NA, format = "fasta", ...)
{
  if (missing(x) || !is(x, "ModStringSet"))
    stop("'x' must be an ModStringSet object", call. = FALSE)
  if (missing(filepath) || !.is_a_string(filepath))
    stop("'filepath' must be a single character value.", call. = FALSE)
  if (!.is_a_string(format))
    stop("'format' must be a single character value.", call. = FALSE)
  if (!.is_a_bool(append))
    stop("'append' must be TRUE or FALSE.", call. = FALSE)
  format <- match.arg(tolower(format), c("fasta", "fastq"))
  compress <- .normarg_compress(compress)
  filepath2 <- path.expand(filepath)
  conFUN <- switch(compress,
                   no = file,
                   gzip = gzfile,
                   bzip2 = bzfile,
                   xz = xzfile)
  con <- conFUN(description = filepath2,
                open = ifelse(append,"ab","wb"),
                encoding = "UTF-8")
  res <- try(switch(format,
                    "fasta" = .write_ModStringSet_to_fasta(x, con, ...),
                    "fastq" = .write_ModStringSet_to_fastq(x, con, ...)),
             silent = TRUE)
  close(con)
  if (is(res, "try-error")) {
    if (!append) {
      if (!file.remove(filepath2))
        warning("cannot remove file '", filepath2, "'")
    }
    stop(res,call. = FALSE)
  }
  invisible(NULL)
}
