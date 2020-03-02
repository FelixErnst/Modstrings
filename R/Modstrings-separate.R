#' @include Modstrings.R
NULL

.from_ModStringSet_to_GRanges <- function(stringset, qualities, codec,
                                          nc.type = "short")
{
  nlength <- length(stringset)
  # we need names for factor generation latter on
  if(is.null(names(stringset))){
    names(stringset) <- seq_along(stringset)
  }
  modInfo <- additionalInfo(codec)
  modInfo <- modInfo[modInfo$abbrev %in% letters(codec),]
  pattern <- letters(codec)[match(modInfo$abbrev, letters(codec))]
  # presearch
  searchstring <- unlist(stringset)
  
  f <- stringi::stri_detect_fixed(as.character(searchstring),pattern)
  # search for modifications
  MIndexList <- lapply(pattern[f],
                       vmatchPattern,
                       stringset)
  names(MIndexList) <- pattern[f]
  f_found <- vapply(MIndexList,
                    function(mi){
                      any(lengths(width(mi)) > 0L)
                    },
                    logical(1))
  if(sum(f_found) != sum(f)){
    stop("Something went wrong.")
  }
  nc_names <- switch(nc.type,
                     "short" = modInfo$short_name[f],
                     "nc" = modInfo$nc[f])
  # construct GRanges object
  ranges <- lapply(MIndexList,unlist)
  f_names <- factor(unlist(lapply(ranges,names)))
  ranges <- unname(unlist(IRanges::IRangesList(ranges)))
  names(ranges) <- f_names
  sns <- mapply(
    function(mi,ncn){
      rep(ncn,sum(lengths(mi)))
    },
    MIndexList,
    nc_names,
    SIMPLIFY = FALSE)
  mcols <- S4Vectors::DataFrame(mod = unlist(sns))
  ans <- unname(GenomicRanges::GRanges(seqnames = names(ranges),
                                       ranges = ranges,
                                       strand = "+",
                                       mcols))
  if(any(!is.na(unlist(qualities)))){
    qualities <- matrix(unlist(qualities),nrow = nlength,byrow = TRUE)
    f_qualities <- matrix(unlist(lapply(coverage(split(ranges,f_names)),
                                        as.logical)),
                          nrow = nlength,
                          byrow = TRUE)
    ans$quality <- t(qualities)[t(f_qualities)]
  }
  ans <- ans[order(match(seqnames(ans),names(stringset)),start(ans))]
  ans
}

#' @rdname separate
#' @aliases combineIntoModstrings separate removeIncompatibleModifications
#'
#' @title Separating and combining a modification information into/from a
#'   \code{XString} and a \code{GRanges} object
#' 
#' @description 
#' With \code{combineIntoModstrings} and \code{separate} the construction and 
#' deconstruction of ModString Objects from an interacive session avoiding 
#' problematic encoding issues. In addition, modification information can be 
#' transfered from/to tabular data with these functions.
#' 
#' \code{combineIntoModstrings} expects \code{seqnames(gr)} or \code{names(gr)} 
#' to match the available \code{names(x)}. Only information with strand 
#' information \code{*} and \code{+} are used.
#' 
#' \code{separate} when used with a \code{GRanges}/\code{GRangesList} object
#' will return an object of the same type, but with modifications seperated. For
#' example an element with \code{mod = "m1Am"} will be returned as two elements
#' with \code{mod = c("Am","m1A")}. The reverse operation is available via
#' \code{combineModifications()}.
#' 
#' \code{removeIncompatibleModifications} filters incompatible modification from
#' a \code{GRanges} or \code{GRangesList}. \code{incompatibleModifications()} 
#' returns the logical vector used for this operation.
#'
#' @param x For \code{separate}: a \code{ModString}/\code{ModStringSet} or 
#' \code{GRanges}/\code{GRangesList}object
#' 
#' For \code{combineIntoModstrings}: a \code{XString} and a \code{XStringSet} 
#' object.
#' @param gr a GRanges object
#' @param with.qualities \code{TRUE} or \code{FALSE} (default): Should the 
#' values from a \code{score} column of the \code{GRanges} object stored? If
#' set \code{with.qualities = TRUE}, \code{combineIntoModstrings} will try to
#' construct a \code{\link{QualityScaledModStringSet}} object.
#' @param quality.type the type of \code{QualityXStringSet} used, if 
#' \code{with.qualities = TRUE}. Must be on of the following values: 
#' \code{"Phred","Solexa","Illumina"}.
#' @param nc.type the type of nomenclature to be used. Either "short" or "nc".
#' "Short" for m3C would be "m3C", "nc" for m3C would be "3C". (
#' \code{default = "short"})
#' @param stop.on.error For \code{combineIntoModstrings}: \code{TRUE}(default)
#' or \code{FALSE}: Should an error be raised upon encounter of incompatible 
#' positions?
#' @param verbose For \code{combineIntoModstrings}: \code{TRUE} or \code{FALSE}
#' (default): Should verbose information reported on the positions filled with 
#' modifications? This settings is passed onto \code{\link{modifyNucleotides}}.
#' @param ... 
#' \itemize{
#' \item{\code{default.quality}:} {For \code{combineIntoModstrings}: the 
#' \code{default.quality} default value for non-modified positions. (default: 
#' \code{default.quality = 0L})}
#' } 
#'
#' @return for \code{separate} a \code{GRanges} object and for 
#' \code{combineIntoModstrings} a \code{ModString*} object or a 
#' \code{QualityScaledModStringSet}, if \code{with.qualities = TRUE}.
#' 
#' @export
#'
#' @examples
#' library(GenomicRanges)
#' # ModDNAString
#' seq <- ModDNAString(paste(alphabet(ModDNAString()), collapse = ""))
#' seq
#' gr <- separate(seq)
#' gr
#' seq2 <- combineIntoModstrings(as(seq,"DNAString"),gr)
#' seq2
#' seq == seq2
#' # ModRNAString
#' seq <- ModRNAString(paste(alphabet(ModRNAString()), collapse = ""))
#' seq
#' gr <- separate(seq)
#' gr
#' # Separating RNA modifications
#' gr <- gr[1]
#' separate(gr)
#' # ... and combine them again (both operations work only on a subset of
#' # modifications)
#' combineModifications(separate(gr))
#' 
#' # handling incompatible modifications
#' seq <- RNAString("AGCU")
#' gr <- GRanges(c("chr1:1:+","chr1:2:+"),mod="m1A")
#' incompatibleModifications(gr,seq)
#' #
#' removeIncompatibleModifications(gr,seq)
setMethod(
  "separate",
  signature = "ModString",
  function(x, nc.type = c("short","nc"))
  {
    nc.type <- match.arg(nc.type)
    codec <- modscodec(seqtype(x))
    ans <- .from_ModStringSet_to_GRanges(as(x, paste0(seqtype(x), "StringSet")),
                                         rep(NA,length(x)), codec = codec,
                                         nc.type = nc.type)
    ans 
  }
)
#' @rdname separate
#' @export
setMethod(
  "separate",
  signature = "ModStringSet",
  function(x, nc.type = c("short","nc"))
  {
    nc.type <- match.arg(nc.type)
    if(is(x,"QualityScaledModStringSet")){
      qualities <- quality(x)
      qualities <- as(qualities,"IntegerList")
    } else {
      qualities <- lapply(width(x),function(i)rep(NA,i))
    }
    codec <- modscodec(seqtype(x))
    ans <- .from_ModStringSet_to_GRanges(x, qualities, codec = codec,
                                         nc.type = nc.type)
    ans
  }
)

# separate for GRanges ---------------------------------------------------------

.norm_GRanges_for_separate <- function(gr)
{
  if(!("mod" %in% colnames(S4Vectors::mcols(gr)))){
    stop("GRanges object does not contain a 'mod' column.",
         call. = FALSE)
  }
  if(any(unique(width(gr)) != 1)){
    stop("width() of GRanges elements must all be == 1",
         call. = FALSE)
  }
  gr
}

.norm_GRangesList_for_separate <- function(grl)
{
  relist(.norm_GRanges_for_separate(unlist(grl, use.names = FALSE)), grl)
}

.detect_seqtype_in_GRanges <- function(gr){
  mod <- mcols(gr)$mod
  if(all(mod %in% names(modsshortnames("ModRNA")))){
    return("ModRNA")
  }
  if(all(mod %in% names(modsshortnames("ModDNA")))){
    return("ModDNA")
  }
  stop("Couldn't detect the seqtype of value in 'mod' column.")
}

.detect_seqtype_in_GRangesList <- function(grl){
  .detect_seqtype_in_GRanges(unlist(grl, use.names = FALSE))
}

# create a library for converting combined modifications (as names) to seperate 
# modifications (as values)
.get_short_nc_lib <- function(seqtype){
  short_nc <- names(modsnomenclature(seqtype))
  # split into codes and the base
  base <- regmatches(short_nc,regexpr("[AGCUX]$",short_nc))
  code <- regmatches(short_nc,regexpr("^[^AGCUX]+",short_nc))
  code_split <- strsplit(code,"")
  # methylattion are always first
  ribo_meth <- vapply(code_split,"[",character(1),1L) == "0" &
    vapply(code_split,"[",character(1),2L) != "0"
  ribo_meth[is.na(ribo_meth)] <- FALSE
  # assemble result lib
  lib <- as.list(short_nc)
  # add the ribometh variations
  lib[ribo_meth] <- pc(paste0("0",base[ribo_meth]),
                       paste0(substr(code[ribo_meth],2,length(code[ribo_meth])),
                              base[ribo_meth]))
  #
  names(lib) <- short_nc
  lib
}

# try to separate modification annotation. Only possible for RNA modifications
.separate_modifications <- function(gr, seqtype)
{
  nc <- .get_nc_ident(mcols(gr)$mod, seqtype, "short")
  lib <- .get_short_nc_lib(seqtype)
  new_nc <- lib[nc]
  gr <- gr[unlist(Map(rep,seq_along(gr),lengths(new_nc)))]
  f <- match(unlist(new_nc), names(modsnomenclature(seqtype)))
  mcols(gr)$mod <- names(modsshortnames(seqtype))[f]
  gr
}
.separate_modifications_in_GRanges <- function(gr, seqtype){
  .separate_modifications(gr, seqtype)
}
.separate_modifications_in_GRangesList <- function(grl, seqtype){
  names <- names(grl)
  if(is.null(names)){
    names(grl) <- seq_along(grl)
  }
  unlisted_gr <- unlist(grl)
  unlisted_gr <- .separate_modifications(unlisted_gr, seqtype)
  ans <- split(unname(unlisted_gr),
               factor(names(unlisted_gr),unique(names(unlisted_gr))))
  names(ans) <- names
  ans
}

#' @rdname separate
#' @export
setMethod(
  "separate",
  signature = "GRanges",
  function(x)
  {
    x <- .norm_GRanges_for_separate(x)
    seqtype <- .detect_seqtype_in_GRanges(x)
    if(seqtype != "ModRNA"){ # currently only supported for "ModRNA"
      return(x)
    }
    x <- .separate_modifications_in_GRanges(x, seqtype)
    x 
  }
)
#' @rdname separate
#' @export
setMethod(
  "separate",
  signature = "GRangesList",
  function(x)
  {
    x <- .norm_GRangesList_for_separate(x)
    seqtype <- .detect_seqtype_in_GRangesList(x)
    if(seqtype != "ModRNA"){ # separate currently only supported for "ModRNA"
      return(x)
    }
    x <- .separate_modifications_in_GRangesList(x, seqtype)
    x
  }
)

# combine ----------------------------------------------------------------------

# remove any unneccessary information
# warn if minus strand information are present
.norm_GRanges <- function(gr, drop.additional.columns = TRUE)
{
  mcols <- mcols(gr)
  if(drop.additional.columns){
    mcols <- mcols[, colnames(mcols) %in% c("mod","quality"), drop = FALSE]
  }
  gr <- GRanges(seqnames = as.character(seqnames(gr)),
                ranges = ranges(gr),
                strand = strand(gr),
                mcols)
  if(any(as.character(strand(gr)) == "-")){
    warning("Annotations on the minus strand will not be used.",
            call. = FALSE)
  }
  gr[as.character(strand(gr)) %in% c("*","+")]
}

.norm_GRangesList <- function(gr, drop.additional.columns = TRUE)
{
  ans <- unlist(gr, use.names = FALSE)
  ans <- .norm_GRanges(ans, drop.additional.columns)
  split(ans, seqnames(ans))
}


# trying to combine modifications if possible ----------------------------------

# return the nc independent whether shortName or nomenclature was used
.get_nc_ident <- function(mod, seqtype, nc_type){
  if(nc_type == "short"){
    modNames <- names(modsshortnames(seqtype))
    f <- match(mod,modNames)
    return(names(modsnomenclature(seqtype))[f])
  }
  if(nc_type == "nc"){
    return(mod)
  }
  stop("Something went wrong.", call. = FALSE)
}

.merge_non_unique_mcols_do <- function(col){
  if(is(col,"List")){
    col <- unlist(col, use.names = FALSE)
  }
  return(relist(col,IRanges::PartitioningByEnd(length(col))))
}

.merge_non_unique_mcols <- function(mcols){
  unique_mcols <- lapply(mcols,unlist)
  unique_mcols <- lapply(unique_mcols,unique)
  f_unique <- lengths(unique_mcols) == 1L
  ans <- mcols[1L,,drop=FALSE]
  if(any(!f_unique)){
    ans[,!f_unique] <- lapply(mcols[,!f_unique,drop=FALSE],
                              .merge_non_unique_mcols_do)
  }
  ans
}

# assembles a new modification from nomenclature and checks if the new type
# is valid
.combine_to_new_nc_ident <- function(gr, seqtype)
{
  nc_type <- .get_nc_type(gr$mod,seqtype)
  nc <- .get_nc_ident(gr$mod,seqtype,nc_type)
  if(is.null(nc)){
    return(gr)
  }
  if(length(unique(nc)) == 1L){
    return(unique(gr))
  }
  # this enforces a hierachy in the nomenclature the smallest numeric value is
  # alsways put first
  nc <- sort(nc)
  tmp <- IRanges::CharacterList(strsplit(nc,""))
  f <- tmp %in% (seq_len(10)-1L)
  base <- unique(unlist(tmp[!f]))
  if(length(base) != 1L){
    stop("Multiple modifications found for position '",unique(start(gr)),"'. ",
         "They could not be combined into a single modification since the base",
         " nucleotides do not match.",
         call. = FALSE)
  }
  newnc <- vapply(tmp[f],
                  paste0,
                  character(1),
                  collapse = "")
  # order the nc base on first digit
  newnc <- newnc[order(vapply(tmp[f],"[",character(1),1L))]
  newnc <- paste0(paste0(newnc,collapse = ""),base)
  newnc <- newnc[newnc %in% names(modsnomenclature(seqtype))]
  if(length(newnc) == 0L){
    stop("Multiple modifications found for position '",unique(start(gr)),"'.",
         "They could not be combined into a single modification.",
         call. = FALSE)
  }
  if(nc_type == "short"){
    f <- match(newnc,names(modsnomenclature(seqtype)))
    newnc <- names(modsshortnames(seqtype))[f]
  }
  mcols <- mcols(gr)
  mcols$mod <- newnc
  mcols <- .merge_non_unique_mcols(mcols)
  GRanges(seqnames = unique(as.character(seqnames(gr))),
          ranges = unique(ranges(gr)),
          strand = unique(as.character(strand(gr))),
          mcols)
}

# sort in new modification annotation
.combine_modifications <- function(gr, seqtype)
{
  f_dup <- duplicated(gr)
  if(any(f_dup)){
    overlappingPos <- gr[f_dup]
    overlappingPos <- unique(overlappingPos)
    f_gr <- gr %in% overlappingPos
    newgr <- lapply(split(gr[f_gr],factor(as.character(gr[f_gr]),
                                          unique(as.character(gr[f_gr])))),
                    .combine_to_new_nc_ident,
                    seqtype)
    newgr <- unlist(GenomicRanges::GRangesList(unname(newgr)))
    gr <- c(gr[!f_gr],newgr)
    gr <- gr[order(gr)]
  }
  gr
}
.combine_modifications_in_GRanges <- function(gr,seqtype){
  .combine_modifications(gr,seqtype)
}
.combine_modifications_in_GRangesList <- function(gr, seqtype){
  unlisted_gr <- unlist(gr)
  unlisted_gr <- .combine_modifications(unlisted_gr, seqtype)
  split(unname(unlisted_gr),
        factor(names(unlisted_gr),unique(names(unlisted_gr))))
}

# normalize GRanges inputs and check compatibility for combineIntoModstrings
.norm_GRanges_for_combine <- function(x, gr, drop.additional.columns = TRUE,
                                      do.combine = FALSE)
{
  gr <- .norm_GRanges(gr, drop.additional.columns)
  if(is(x,"XString")){
    if(any(max(end(gr)) > length(x))){
      stop("GRanges object contains coordinates out of bounds for the ",
           "XStringSet object.",
           call. = FALSE)
    }
  } else {
    x <- do.call(paste0(x,"String"),list("N"))
  }
  if(any(unique(width(gr)) != 1)){
    stop("width() of GRanges elements must all be == 1",
         call. = FALSE)
  }
  if(!("mod" %in% colnames(S4Vectors::mcols(gr)))){
    stop("GRanges object does not contain a 'mod' column.",
         call. = FALSE)
  }
  # check if modifications are compatible with type of input string
  seqtype <- paste0("Mod",gsub("Mod","",seqtype(x)))
  .get_nc_type(gr$mod, seqtype)
  # this can also be done with checking the findOverlaps(gr) length. This
  # should however be faster
  f <- duplicated(start(gr))
  if(any(f) & do.combine){
    f <- which(duplicated(start(gr)))
    msg <- paste0("Multiple modifications found for position '", start(gr)[f],
                  "'.")
    if(!is.null(mcols(gr)$quality)){
      stop(msg, " Remove quality metadata column or remove overlapping",
           " modifications.", call. = FALSE)
    }
    gr <- .combine_modifications_in_GRanges(gr,seqtype)
    if(any(duplicated(start(gr)))){
      stop(msg, " They could not be combined.", call. = FALSE)
    }
  }
  gr
}

.norm_GRangesList_for_combine <- function(x, gr, drop.additional.columns = TRUE,
                                          do.combine = FALSE)
{
  gr <- .norm_GRangesList(gr, drop.additional.columns)
  if(length(gr) == 0L){
    stop("'gr' is empty (length() == 0L)", call. = FALSE)
  }
  if(!("mod" %in% colnames(S4Vectors::mcols(gr@unlistData)))){
    stop("GRanges object does not contain a 'mod' column.",
         call. = FALSE)
  }
  #
  if(is(x,"XStringSet")){
    seqnames <- as(seqnames(gr),"CharacterList")
    if(is.null(names(x)) ||
       !all(unique(unlist(seqnames)) %in% names(x))){
      stop("Names of XStringSet object do not match the seqnames or names in the",
           " elements of the GRangesList object.",
           call. = FALSE)
    }
    f <- names(x) %in% unlist(unique(seqnames))
    m <- match(names(x), unlist(unique(seqnames)))
    m <- m[!is.na(m)]
    if(any(max(end(gr))[m] > width(x[f]))){
      stop("Elements of the GRangesList object contain coordinates out of bounds",
           " for the XStringSet object.",
           call. = FALSE)
    }
    if(any(unique(unlist(width(gr))) != 1)){
      stop("width() of GRangesList elements must all be == 1",
           call. = FALSE)
    }
  } else {
    x <- do.call(paste0(x,"StringSet"),list("N"))
  }
  # check if modifications are compatible with type of input string
  seqtype <- paste0("Mod",gsub("Mod","",seqtype(x)))
  .get_nc_type(unlist(S4Vectors::mcols(gr@unlistData)$mod), seqtype)
  #
  starts <- start(gr)
  starts <- any(duplicated(starts))
  if(any(starts) & do.combine){
    msg <- paste0("Multiple modifications found for position in element '",
                  which(starts),"'.")
    qualityCol <- "quality" %in% colnames(unlist(S4Vectors::mcols(gr,level="within")))
    if(any(qualityCol)){
      stop(msg," Remove quality metadata column or remove overlapping",
           " modifications.",
           call. = FALSE)
    }
    gr <- .combine_modifications_in_GRangesList(gr, seqtype)
    starts <- start(gr)
    starts <- vapply(starts, function(z){any(duplicated(z))},logical(1))
    if(any(starts)){
      stop(msg, " They could not be combined.", call. = FALSE)
    }
  }
  gr
}

# check which annotation types is it. If nothing is found the seqtype is not
# compatible with the mod type
.get_nc_type <- function(mod,seqtype)
{
  sn <- names(modsshortnames(seqtype))
  nc <- names(modsnomenclature(seqtype))
  if(all(mod %in% sn)){
    return("short")
  }
  if(all(mod %in% nc)){
    return("nc")
  }
  if(any(mod %in% sn) || any(mod %in% nc)){
    stop("Nomenclature type of the modifications must be unique. The 'GRanges'",
         "object contains both types.", call. = FALSE)
  }
  stop("Nomenclature types of modifications are not compatible with '",
       seqtype,"String'.",
       call. = FALSE)
}

# convert the position information into a logical list or matrix
.pos_to_logical_list <- function(x, at)
{
  width <- width(x)
  list <- lapply(width,function(w){rep(FALSE,w)})
  if(!is.null(names(x)) && !is.null(names(at))){
    names(list) <- names(x)
    f <- names(x) %in% names(at)
    list <- list[f]
  }
  if(length(list) != length(at)){
    stop("Length of 'x' and 'gr' does not match.", call. = FALSE)
  }
  list <- mapply(
    function(l,j){
      l[j] <- TRUE
      l
    },
    list,
    at,
    SIMPLIFY = FALSE)
  list
}

.pos_to_logical_matrix <- function(x, at)
{
  width <- width(x)
  m <- matrix(rep(FALSE,sum(width)),length(x))
  if(!is.null(names(x)) && !is.null(names(at))){
    f <- names(x) %in% names(at)
    m <- m[f,,drop=FALSE]
  }
  if(nrow(m) != length(at)){
    stop("Length of 'x' and 'gr' does not match.", call. = FALSE)
  }
  m <- mapply(
    function(i,j){
      n <- m[i,,drop = FALSE]
      n[,j] <- TRUE
      n
    },
    seq_len(nrow(m)),
    at,
    SIMPLIFY = FALSE)
  m <- matrix(unlist(m),length(m),byrow = TRUE)
  m
}

.get_XStringQuality_from_GRanges <- function(string, gr, quality.type,
                                             default.quality)
{
  if(is(string,"ModString")){
    width  <- length(string)
    pos <- list(start(gr[[1]]))
  } else if(is(string,"ModStringSet")){
    width  <- width(string)
    pos <- as.list(start(gr))
  }
  qualities <- lapply(width, function(i) rep(default.quality,i))
  modQualities <- lapply(gr, function(g){S4Vectors::mcols(g)$quality})
  f <- names(string) %in% names(gr)
  qualities[f] <- mapply(
    function(q,p,mq){
      q[p] <- mq
      q
    },
    qualities[f],
    pos,
    modQualities,
    SIMPLIFY = FALSE)
  ans <- do.call(paste0(quality.type,"Quality"),
                 list(as(qualities,"IntegerList")))
  ans
}

# construct a QualityScaleStringSet from information in the GRanges object
.convert_to_QualityScaledStringSet <- function(string, gr, quality.type, ...)
{
  qualityCol <- vapply(gr,
                       function(g){
                         "quality" %in% colnames(S4Vectors::mcols(g))
                       },
                       logical(1))
  if(!all(qualityCol)){
    stop("If with.qualities is set to FALSE, all GRanges objects have to ",
         "contain a 'quality' column.",
         call. = FALSE)
  }
  args <- list(...)
  if(!is.null(args[["default.quality"]])){
    default.quality <- as.integer(args[["default.quality"]])
  } else {
    default.quality <- 0L
  }
  qualities <- .get_XStringQuality_from_GRanges(string, gr, quality.type,
                                                default.quality)
  do.call(paste0("QualityScaled",seqtype(string),"StringSet"),
          list(string, qualities))
}

#' @rdname separate
#' @export
setMethod(
  "combineIntoModstrings",
  signature = c(x = "XString", gr = "GRanges"),
  function(x, gr, with.qualities = FALSE, quality.type = "Phred",
           stop.on.error = TRUE, verbose = FALSE, ...)
  {
    .check_verbose(verbose)
    .check_stop.on.error(stop.on.error)
    x <- as(x,
            paste0("Mod",
                   seqtype(x),
                   "String"))
    gr <- .norm_GRanges_for_combine(x, gr, do.combine = TRUE)
    nc.type <- .get_nc_type(gr$mod,seqtype(x))
    at <- .pos_to_logical_matrix(as(x, paste0(seqtype(x), "StringSet")),
                                 list(start(gr)))
    ans_seq <- modifyNucleotides(x, as.vector(at), S4Vectors::mcols(gr)$mod,
                                 nc.type = nc.type, 
                                 stop.on.error = stop.on.error,
                                 verbose = verbose)
    if(!with.qualities){
      return(ans_seq)
    }
    quality.type <- match.arg(quality.type,c("Phred", "Solexa", "Illumina"))
    .convert_to_QualityScaledStringSet(ans_seq, gr, quality.type, ...)
  }
)

#' @rdname separate
#' @export
setMethod(
  "combineIntoModstrings",
  signature = c(x = "XStringSet", gr = "GRangesList"),
  function(x, gr, with.qualities = FALSE, quality.type = "Phred",
           stop.on.error = TRUE, verbose = FALSE, ...)
  {
    .check_verbose(verbose)
    .check_stop.on.error(stop.on.error)
    if(!assertive::is_a_bool(with.qualities)){
      stop("with.qualities has to be TRUE or FALSE.")
    }
    quality.type <- match.arg(quality.type, c("Phred", "Solexa", "Illumina"))
    gr <- .norm_GRangesList_for_combine(x, gr, do.combine = TRUE)
    m <- match(names(x),names(gr))
    f <- !is.na(m)
    m <- m[f]
    if (!S4Vectors::isConstant(width(x))){
      at <- .pos_to_logical_list(as(x, paste0(seqtype(x), "StringSet"))[f],
                                 start(gr)[m])
    } else {
      at <- .pos_to_logical_matrix(as(x, paste0(seqtype(x), "StringSet"))[f],
                                   start(gr)[m])
    }
    mod <- S4Vectors::mcols(gr[m], level="within")[,"mod"]
    ans_seq <- modifyNucleotides(x[f], at, mod,
                                 stop.on.error = stop.on.error,
                                 verbose = verbose)
    ans_seq <- c(ans_seq,x[!f])
    ans_seq <- ans_seq[match(names(x),names(ans_seq))]
    if(!with.qualities){
      return(ans_seq)
    }
    .convert_to_QualityScaledStringSet(ans_seq, gr, quality.type, ...)
  }
)

#' @rdname separate
#' @export
setMethod(
  "combineIntoModstrings",
  signature = c(x = "XStringSet", gr = "GRanges"),
  function(x, gr, with.qualities = FALSE, quality.type = "Phred",
           stop.on.error = TRUE, verbose = FALSE, ...)
  {
    .check_verbose(verbose)
    .check_stop.on.error(stop.on.error)
    if(!assertive::is_a_bool(with.qualities)){
      stop("with.qualities has to be TRUE or FALSE.")
    }
    quality.type <- match.arg(quality.type, c("Phred", "Solexa", "Illumina"))
    # If names are set use these. Otherwise use seqnames
    if(is.null(names(gr))){
      gr <- split(gr, seqnames(gr))
      gr <- gr[lengths(gr) != 0L]
    } else {
      gr <- split(unname(gr), names(gr))
    }
    gr <- .norm_GRangesList_for_combine(x, gr, do.combine = TRUE)
    m <- match(names(x),names(gr))
    f <- !is.na(m)
    m <- m[f]
    if (!S4Vectors::isConstant(width(x))){
      at <- .pos_to_logical_list(as(x, paste0(seqtype(x), "StringSet"))[f],
                                 start(gr)[m])
    } else {
      at <- .pos_to_logical_matrix(as(x, paste0(seqtype(x), "StringSet"))[f],
                                   start(gr)[m])
    }
    mod <- S4Vectors::mcols(gr[m], level="within")[,"mod"]
    ans_seq <- modifyNucleotides(x[f], at, mod,
                                 stop.on.error = stop.on.error,
                                 verbose = verbose)
    ans_seq <- c(ans_seq,x[!f])
    ans_seq <- ans_seq[match(names(x),names(ans_seq))]
    if(!with.qualities){
      return(ans_seq)
    }
    .convert_to_QualityScaledStringSet(ans_seq, gr, quality.type, ...)
  }
)

# combineModifications ---------------------------------------------------------

#' @rdname separate
#' @export
setMethod(
  "combineModifications",
  signature = c(gr = "GRanges"),
  function(gr)
  {
    seqtype <- .detect_seqtype_in_GRanges(gr)
    gr <- .norm_GRanges_for_combine(seqtype, gr, drop.additional.columns = FALSE,
                                    do.combine = TRUE)
    gr
  }
)

#' @rdname separate
#' @export
setMethod(
  "combineModifications",
  signature = c(gr = "GRangesList"),
  function(gr)
  {
    seqtype <- .detect_seqtype_in_GRanges(gr)
    gr <- .norm_GRangesList_for_combine(seqtype, gr, drop.additional.columns = FALSE,
                                        do.combine = TRUE)
    gr
  }
)


# incompatibleModifications ---------------------------------------------------- 

.incompatbile_modifications <- function(x, at, mod){
  at <- .norm_replace_pos_ModString(x,at)
  assertive::assert_all_are_non_empty_character(as.character(unlist(mod)))
  if(length(at) != length(mod)){
    stop("lengths of 'at' and 'mod' need to be equal.",
         call. = FALSE)
  }
  # check if originating base matches the modification
  seqtype <- paste0("Mod",gsub("Mod","",seqtype(x)))
  modValues <- .norm_seqtype_modtype(unlist(mod), seqtype, "short", class(x))
  codec <- modscodec(seqtype)
  f <- values(codec)[match(modValues, values(codec))]
  current_letter <- lapply(at,
                           function(i){
                             subseq(x,i,i)
                           })
  current_letter <- vapply(current_letter,as.character,character(1))
  if(is(x,"ModString")){
    class <- paste0(class(x),"Set")
    current_letter <- as(do.call(class, list(current_letter)),
                         gsub("Mod","",class))
    current_letter <- as.character(current_letter)
  }
  originating_base <- lapply(relist(f,mod),
                             function(i){originatingBase(codec)[i]})
  mismatch <- IRanges::CharacterList(originating_base) != current_letter
  unlist(mismatch, use.names = FALSE)
}

#' @rdname separate
#' @export
setMethod(
  "incompatibleModifications",
  signature = c(gr = "GRanges", x = "XString"),
  function(gr, x)
  {
    gr <- .norm_GRanges_for_combine(x, gr, drop.additional.columns = FALSE,
                                    do.combine = FALSE)
    at <- .pos_to_logical_matrix(as(x, paste0(seqtype(x), "StringSet")),
                                 list(start(gr)))
    mod <- S4Vectors::mcols(gr)$mod
    mod <- split(mod, start(gr))
    mismatch <- .incompatbile_modifications(x, as.vector(at), mod)
    mismatch
  }
)


#' @rdname separate
#' @export
setMethod(
  "incompatibleModifications",
  signature = c(gr = "GRanges", x = "XStringSet"),
  function(gr, x)
  {
    if(is.null(names(gr))){
      gr <- split(gr, seqnames(gr))
      gr <- gr[lengths(gr) != 0L]
    } else {
      gr <- split(unname(gr), names(gr))
    }
    unlist(incompatibleModifications(gr, x), use.names = FALSE)
  }
)

#' @rdname separate
#' @export
setMethod(
  "incompatibleModifications",
  signature = c(gr = "GRangesList", x = "XStringSet"),
  function(gr, x)
  {
    gr <- .norm_GRangesList_for_combine(x, gr, drop.additional.columns = FALSE,
                                        do.combine = FALSE)
    m <- match(names(x),names(gr))
    f <- !is.na(m)
    m <- m[f]
    at <- .pos_to_logical_list(as(x, paste0(seqtype(x), "StringSet"))[f],
                               start(gr)[m])
    mod <- S4Vectors::mcols(gr[m], level="within")[,"mod"]
    mismatch <- Map(.incompatbile_modifications,
                    x[f],
                    at,
                    mod)
    mismatch <- relist(unlist(mismatch, use.names = FALSE), gr)
    mismatch
  }
)

# removeIncompatibleModifications ----------------------------------------------

#' @rdname separate
#' @export
setMethod(
  "removeIncompatibleModifications",
  signature = c(gr = "GRanges", x = "XString"),
  function(gr, x)
  {
    gr <- .norm_GRanges_for_combine(x,gr, drop.additional.columns = FALSE,
                                    do.combine = FALSE)
    at <- .pos_to_logical_matrix(as(x, paste0(seqtype(x), "StringSet")),
                                 list(start(gr)))
    mod <- S4Vectors::mcols(gr)$mod
    mod <- split(mod, start(gr))
    mismatch <- .incompatbile_modifications(x, as.vector(at), mod)
    gr[!mismatch]
  }
)

#' @rdname separate
#' @export
setMethod(
  "removeIncompatibleModifications",
  signature = c(gr = "GRanges", x = "XStringSet"),
  function(gr, x)
  {
    if(is.null(names(gr))){
      gr <- split(gr, seqnames(gr))
      gr <- gr[lengths(gr) != 0L]
    } else {
      gr <- split(unname(gr), names(gr))
    }
    unlist(removeIncompatibleModifications(gr, x), use.names = FALSE)
  }
)

#' @rdname separate
#' @export
setMethod(
  "removeIncompatibleModifications",
  signature = c(gr = "GRangesList", x = "XStringSet"),
  function(gr, x)
  {
    gr <- .norm_GRangesList_for_combine(x, gr, drop.additional.columns = FALSE,
                                        do.combine = FALSE)
    m <- match(names(x),names(gr))
    f <- !is.na(m)
    m <- m[f]
    at <- .pos_to_logical_list(as(x, paste0(seqtype(x), "StringSet"))[f],
                               start(gr)[m])
    mod <- S4Vectors::mcols(gr[m], level="within")[,"mod"]
    mod <- lapply(mod,split,start(gr))
    mismatch <- Map(.incompatbile_modifications,
                    x[f],
                    at,
                    mod)
    mismatch <- relist(unlist(mismatch, use.names = FALSE), gr)
    gr <- gr[!mismatch]
    gr[lengths(gr) != 0L]
  }
)
