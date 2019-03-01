#' @include Modstrings.R
#' @include Modstrings-ModStringSet.R
NULL

#' @name ModStringViews
#' @aliases ModStringViews show,ModStringViews-method 
#' ==,ModStringViews,ModStringViews-method ==,ModStringViews,XString-method 
#' ==,XStringViews,ModString-method ModStringSet,ModStringViews-method
#' 
#' @title The ModStringViews class extending the XStringViews class
#' 
#' @description 
#' As the \code{\link[Biostrings:XStringViews-class]{XStringViews}} the 
#' \code{ModStringViews} is the basic container for storing a set of views on 
#' the same sequence (this time a \code{ModString} object).
#' 
#' @details 
#' For the details have a look at the 
#' \code{\link[Biostrings:XStringViews-class]{XStringViews}} class.
#' 
#' @param subject,start,end,width,names
#' See \code{\link[Biostrings:XStringViews-class]{XStringViews}}.
#' @param x a \code{ModStringViews} object.
#' @param row.names passed on to \code{\link[base:as.data.frame]{as.data.frame}}
#' @param optional passed on to \code{\link[base:as.data.frame]{as.data.frame}}
#' @param ... passed on to \code{\link[base:as.data.frame]{as.data.frame}}
#' 
#' @return a \code{ModStringViews} object.
#' 
#' @examples
#' seq <- ModDNAString("AGC6AGC6")
#' v <- Views(seq, start = 3:1, end = 6:8)
NULL

# derived from Biostrings/R/XStringViews-class.R -------------------------------
setClass("ModStringViews",
         contains="XStringViews",
         representation(
           subject = "ModString"
         )
)

# derived from Biostrings/R/XStringViews-class.R -------------------------------
# Constructor

unsafe.newModStringViews <- function(subject, start, width){
  new2("ModStringViews",
       subject = subject,
       ranges = IRanges::IRanges(start = start, width = width),
       check = FALSE)
}
#' @rdname ModStringViews
#' @export
setMethod(
  "Views", "ModString",
  function(subject, start = NULL, end = NULL, width = NULL, 
           names = NULL)
  {
    .new_Views(subject, start = start, end = end, width = width, names = names,
               Class = "ModStringViews")
  }
)
#' @rdname ModStringViews
#' @export
setMethod(
  "Views", "character",
  function(subject, start = NULL, end = NULL, width = NULL, 
           names = NULL)
  {
    xsubject <- ModString(NULL, subject)
    Views(xsubject, start = start, end = end, width = width, names = names)
  }
)

# derived from Biostrings/R/XStringViews-class.R -------------------------------
# Coercion

### We need this so that B/DNA/RNA/AAStringSet() used below work on an
### XStringViews object.

setMethod(
  "ModStringSet",
  "ModStringViews",
  function(seqtype, x, start = NA, end = NA, width = NA, use.names = TRUE)
  {
    y <- .fromXStringViewsToStringSet(x, out.of.limits = "warning",
                                      use.names = use.names)
    
    ModStringSet(seqtype, y, start = start, end = end, width = width,
                 use.names = TRUE)
  }
)

#' @export
setAs("ModStringViews", "ModStringSet",
      function(from) .fromXStringViewsToStringSet(from,
                                                  out.of.limits = "warning",
                                                  use.names = TRUE))
#' @export
setAs("ModStringViews", "ModDNAStringSet", function(from) ModDNAStringSet(from))
#' @export
setAs("ModStringViews", "ModRNAStringSet", function(from) ModRNAStringSet(from))
#' @export
setAs("ModStringSet", "Views", .XStringSetAsViews)
#' @export
setAs("ModStringSet", "ModStringViews", .XStringSetAsViews)
#' @rdname ModStringViews
#' @export
setMethod("as.data.frame", "ModStringViews",
          function (x, row.names = NULL, optional = FALSE, ...){
            as.data.frame(as(x, "ModStringSet"), row.names, optional, ...)
          })
#' @export
setAs("ModStringViews", "ModDNAStringSet", function(from) ModDNAStringSet(from))
#' @export
setAs("ModStringViews", "ModRNAStringSet", function(from) ModRNAStringSet(from))


# derived from Biostrings/R/XStringViews-class.R -------------------------------
# show

# These functions need to be here to access the modified functions of
# - XString.read ==> ModString.read

### nchar(XStringViews.get_view(x, start, end)) is always end-start+1
ModStringViews.get_view <- function(x, start, end)
{
  lx <- length(x)
  if (end < 1 || start > lx)
    return(format("", width = end-start+1))
  Lmargin <- ""
  if (start < 1) {
    Lmargin <- format("", width = 1-start)
    start <- 1
  }
  Rmargin <- ""
  if (end > lx) {
    Rmargin <- format("", width = end-lx)
    end <- lx
  }
  paste0(Lmargin, ModString.read(x, start, end), Rmargin)
}

ModStringViews.get_snippet <- function(x, start, end, snippetWidth)
{
  if (snippetWidth < 7)
    snippetWidth <- 7
  width <- end - start + 1
  if (width <= snippetWidth) {
    ModStringViews.get_view(x, start, end)
  } else {
    w1 <- (snippetWidth - 2) %/% 2
    w2 <- (snippetWidth - 3) %/% 2
    paste0(ModStringViews.get_view(x, start, start+w1-1),"...",
           ModStringViews.get_view(x, end-w2+1, end))
  }
}

ModStringViews.show_vframe_header <- function(iW, startW, endW, widthW)
{
  cat(format("", width = iW+1),
      format("start", width = startW, justify = "right"), " ",
      format("end", width = endW, justify = "right"), " ",
      format("width", width = widthW, justify = "right"), "\n",
      sep = "")
}

ModStringViews.show_vframe_line <- function(x, i, iW, startW, endW, widthW)
{
  start <- start(x)[i]
  end <- end(x)[i]
  width <- end - start + 1
  snippetWidth <- getOption("width") - 6 - iW - startW - endW - widthW
  cat(format(paste0("[", i,"]"), width = iW, justify = "right"), " ",
      format(start, width = startW, justify="right"), " ",
      format(end, width = endW, justify="right"), " ",
      format(width, width = widthW, justify="right"), " ",
      "[", ModStringViews.get_snippet(subject(x), start, end, snippetWidth),
      "]\n",
      sep = "")
}

ModStringViews.show_vframe <- function(x)
{
  nhead <- get_showHeadLines()
  ntail <- get_showTailLines()
  cat("\nviews:")
  lx <- length(x)
  if (lx == 0)
    cat(" NONE\n")
  else {
    cat("\n")
    iW <- nchar(as.character(lx)) + 2 # 2 for the brackets
    startMax <- max(start(x))
    startW <- max(nchar(startMax), nchar("start"))
    endMax <- max(end(x))
    endW <- max(nchar(endMax), nchar("end"))
    widthMax <- max(width(x))
    widthW <- max(nchar(widthMax), nchar("width"))
    ModStringViews.show_vframe_header(iW,
                                    startW,
                                    endW,
                                    widthW)
    if (lx <= nhead + ntail +1) {
      for (i in seq_len(lx))
        ModStringViews.show_vframe_line(x,
                                        i,
                                        iW,
                                        startW,
                                        endW,
                                        widthW)
    } else {
      if (nhead > 0)
        for (i in seq_len(nhead))
          ModStringViews.show_vframe_line(x,
                                          i,
                                          iW,
                                          startW,
                                          endW,
                                          widthW)
      cat(format("...", width = iW, justify = "right"),
          " ",
          format("...", width = startW, justify = "right"),
          " ",
          format("...", width = endW, justify = "right"),
          " ",
          format("...", width = widthW, justify = "right"),
          " ...\n", sep = "")
      if (ntail > 0)
        for (i in (lx-ntail+1L):lx)
          ModStringViews.show_vframe_line(x, i, iW, startW, endW,  widthW)
    }
  }
}

setMethod("show", "ModStringViews",
  function(object)
  {
    subject <- subject(object)
    lsub <- length(subject)
    cat("  Views on a ", lsub, "-letter ", class(subject), " subject", sep = "")
    cat("\nsubject:", .toSeqSnippet(subject, getOption("width") - 9))
    ModStringViews.show_vframe(object)
  }
)


# derived from Biostrings/R/XStringViews-class.R -------------------------------
# Comparison

# These functions need to be here to access the modified functions of
# - comparable_seqtypes

#' @export
setMethod("==", signature(e1 = "ModStringViews", e2 = "ModStringViews"),
          function(e1, e2){
            if (!comparable_seqtypes(seqtype(e1), seqtype(e2))) {
              class1 <- class(subject(e1))
              class2 <- class(subject(e2))
              stop("comparison between XStringViews objects with subjects of ",
                   "class \"", class1, "\" and \"", class2, "\" ",
                   "is not supported")
            }
            .XStringViews.equal(e1, e2)
          }
)

.compare_ModString_Views <- function(e1,
                                     e2){
  if (!comparable_seqtypes(seqtype(e1), seqtype(e2))) {
    class1 <- class(subject(e1))
    class2 <- class(e2)
    stop("comparison between an XStringViews object with a subject of ",
         "class \"", class1, "\" and a \"", class2, "\" instance ",
         "is not supported")
  }
  .XStringViews.equal(e1, as(e2, "Views"))
}
#' @export
setMethod("==", signature(e1 = "ModStringViews", e2 = "XString"), 
          function(e1, e2) .compare_ModString_Views(e1, e2)
)
#' @export
setMethod("==", signature(e1 = "XStringViews", e2 = "ModString"), 
          function(e1, e2) .compare_ModString_Views(e1, e2)
)
