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
#' @title Separating and combining a ModString object into/from a XString and a 
#' GRanges object
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
#' \code{removeIncompatibleModifications} filters incompatible modification from
#' a \code{GRanges} or \code{GRangesList}.
#'
#' @param x For \code{separate}: a \code{ModString} or \code{ModStringSet}
#' object
#' 
#' For \code{combineIntoModstrings}: a \code{XStringSet} and a \code{GRanges} 
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
#' 
#'
#' @return for \code{separate} a \code{GRanges} object and for 
#' \code{combineIntoModstrings} a \code{ModString*} object or a 
#' \code{QualityScaledModStringSet}, if \code{with.qualities = TRUE}.
#' 
#' @export
#'
#' @examples
#' seq <- ModDNAString(paste(alphabet(ModDNAString()), collapse = ""))
#' gr <- separate(seq)
#' seq <- combineIntoModstrings(as(seq,"DNAString"),gr)
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
  newnc <- newnc[order(unlist(tmp[f][,1]))]
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
  GRanges(seqnames = unique(as.character(seqnames(gr))),
          ranges = unique(ranges(gr)),
          strand = unique(as.character(strand(gr))),
          mod = newnc)
}

# sort in new modification annotation
.combine_modifications <- function(gr, seqtype)
{
  overlappingPos <- start(gr)[duplicated(start(gr))]
  if(length(overlappingPos) > 0L){
    newgr <- lapply(overlappingPos,
                    function(i){
                      .combine_to_new_nc_ident(gr[start(gr) == i,],seqtype)
                    })
    gr <- gr[!(start(gr) %in% overlappingPos)]
    gr <- unlist(GRangesList(c(list(gr),newgr)))
    gr <- gr[order(start(gr))]
  }
  gr
}
.combine_modifications_in_GRanges <- function(gr,seqtype){
  .combine_modifications(gr,seqtype)
}
.combine_modifications_in_GRangesList <- function(gr,seqtype){
  gr <- GRangesList(lapply(gr,
                           function(g){
                             .combine_modifications(g,seqtype)
                           }))
  gr
}

# normalize GRanges inputs and check compatibility for combineIntoModstrings
.norm_GRanges_for_combine <- function(x,gr)
{
  gr <- .norm_GRanges(gr)
  if(any(max(end(gr)) > length(x))){
    stop("GRanges object contains coordinates out of bounds for the ",
         "XStringSet object.",
         call. = FALSE)
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
  seqtype <- seqtype(x)
  if(!stringr::str_detect(seqtype,"Mod")){
    seqtype <- paste0("Mod",seqtype)
  }
  .get_nc_type(gr$mod,seqtype)
  # this can also be done with checking the findOverlaps(gr) length. This
  # should however be faster
  f <- duplicated(start(gr))
  if(any(f)){
    if(!is.null(mcols(gr)$quality)){
      stop("Multiple modifications found for position '",start(gr)[f],"'.",
           call. = FALSE)
    }
    gr <- .combine_modifications_in_GRanges(gr,seqtype)
    if(any(duplicated(start(gr)))){
      f <- which(duplicated(start(gr)))
      stop("Multiple modifications found for position '",start(gr)[f],"'.",
           call. = FALSE)
    }
  }
  gr
}

.norm_GRangesList_for_combine <- function(x, gr, drop.additional.columns = TRUE)
{
  gr <- .norm_GRangesList(gr, drop.additional.columns)
  if(length(gr) == 0L){
    stop("'gr' is empty (length() == 0L)", call. = FALSE)
  }
  if(!("mod" %in% colnames(S4Vectors::mcols(gr@unlistData)))){
    stop("GRanges object does not contain a 'mod' column.",
         call. = FALSE)
  }
  # check if modifications are compatible with type of input string
  seqtype <- seqtype(x)
  if(!stringr::str_detect(seqtype,"Mod")){
    seqtype <- paste0("Mod",seqtype)
  }
  .get_nc_type(unlist(S4Vectors::mcols(gr@unlistData)$mod),seqtype)
  #
  seqnames <- as(seqnames(gr),"CharacterList")
  if(is.null(names(x)) ||
     !all(unique(unlist(seqnames)) %in% names(x))){
    stop("Names of XStringSet object do not match the seqnames or names in the",
         " elements of the GRangesList object.",
         call. = FALSE)
  }
  f <- names(x) %in% vapply(seqnames, unique, character(1))
  m <- match(names(x),vapply(seqnames, unique, character(1)))
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
  #
  starts <- start(gr)
  starts <- vapply(starts, function(z){any(duplicated(z))},logical(1))
  if(any(starts)){
    qualityCol <- "quality" %in% colnames(unlist(S4Vectors::mcols(gr,level="within")))
    if(any(qualityCol)){
      stop("Multiple modifications found for position in element '",
           which(starts),"'. Remove quality metadata column or remove ",
           "overlapping modifications.",
           call. = FALSE)
    }
    gr <- .combine_modifications_in_GRangesList(gr,seqtype)
    starts <- start(gr)
    starts <- vapply(starts, function(z){any(duplicated(z))},logical(1))
    if(any(starts)){
      stop("Multiple modifications found for position in element '",
           which(starts),"'.",
           call. = FALSE)
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
    gr <- .norm_GRanges_for_combine(x,gr)
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
    gr <- .norm_GRangesList_for_combine(x, gr)
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
    gr <- .norm_GRangesList_for_combine(x, gr)
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


# removeIncompatibleModifications ----------------------------------------------

.remove_incompatbile_modifications <- function(x, at, mod){
  at <- .check_replace_pos_ModString(x,at)
  assertive::assert_all_are_non_empty_character(as.character(mod))
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
  current_letter <-  as(do.call(paste0(seqtype(current_letter[[1L]]),"StringSet"),
                                list(current_letter)),
                        paste0(gsub("Mod","",seqtype(current_letter[[1L]])),"StringSet"))
  current_letter <- as.character(current_letter)
  mismatch <- originatingBase(codec)[f] != current_letter
  mismatch
}

#' @rdname separate
#' @export
setMethod(
  "removeIncompatibleModifications",
  signature = c(gr = "GRanges", x = "XString"),
  function(gr, x)
  {
    x <- as(x, paste0("Mod", seqtype(x), "String"))
    gr <- .norm_GRanges_for_combine(x,gr)
    at <- .pos_to_logical_matrix(as(x, paste0(seqtype(x), "StringSet")),
                                 list(start(gr)))
    gr[.remove_incompatbile_modifications(x, as.vector(at),
                                          S4Vectors::mcols(gr)$mod)]
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
    gr <- .norm_GRangesList_for_combine(x, gr, drop.additional.columns = FALSE)
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
    part <- IRanges::PartitioningByEnd(gr[m])
    mismatch <- .remove_incompatbile_modifications(unlist(x[f]), unlist(at), 
                                                   unlist(mod))
    ans <- gr[!relist(mismatch,part)]
    ans[lengths(ans) != 0L]
  }
)
