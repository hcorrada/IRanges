### =========================================================================
### DataFrameList objects
### -------------------------------------------------------------------------

setClass("DataFrameList", representation("VIRTUAL"),
         prototype = prototype(elementType = "DataFrame"),
         contains = "Sequence")
setClass("SimpleDataFrameList",
         prototype = prototype(elementType = "DataFrame"),
         contains = c("SimpleList", "DataFrameList"))

setClass("SplitDataFrameList", representation("VIRTUAL"),
         prototype = prototype(elementType = "DataFrame"),
         contains = "DataFrameList")
setClass("SimpleSplitDataFrameList",
         prototype = prototype(elementType = "DataFrame"),
         contains = c("SplitDataFrameList", "SimpleDataFrameList"))
setClass("CompressedSplitDataFrameList",
         prototype = prototype(elementType = "DataFrame",
                               unlistData = new("DataFrame")),
         contains = c("SplitDataFrameList", "CompressedList"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("dim", "SimpleDataFrameList",
          function(x) {
            if (length(x) == 0L)
              c(0L, 0L)
            else {
              nrow <- sum(unlist(lapply(as.list(x, use.names = FALSE), nrow)))
              ncol <- ncol(x[[1]])
              as.integer(c(nrow, ncol))
            }
          })
  
setMethod("dim", "CompressedSplitDataFrameList",
          function(x) {
            if (length(x) == 0L)
              c(0L, 0L)
            else
              dim(unlist(x))
          })

### FIXME: make separate rownames, colnames methods, because the
### rownames calculation can be _very_ slow
setMethod("dimnames", "SimpleDataFrameList",
          function(x) {
            if (length(x) == 0L)
              list(character(), character())
            else
              list(unlist(lapply(x, rownames), use.names=FALSE),
                   colnames(x[[1]]))
          })

setMethod("dimnames", "CompressedSplitDataFrameList",
          function(x) {
            if (length(x) == 0L)
              list(character(), character())
            else
              dimnames(unlist(x))
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.SplitDataFrameList <- function(x) {
  if (length(x)) {
    firstNames <- colnames(x[[1]])
    if (!all(sapply(as.list(x, use.names = FALSE),
                    function(df) identical(firstNames, colnames(df)))))
      return("column counts or names differ across elements")
  }
  NULL
}

setValidity2("SplitDataFrameList", .valid.SplitDataFrameList)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

DataFrameList <- function(...)
{
  SimpleList("SimpleDataFrameList", list(...))
}

SplitDataFrameList <- function(..., compress = FALSE)
{
  if (compress)
    CompressedList("CompressedSplitDataFrameList", list(...))
  else
    SimpleList("SimpleSplitDataFrameList", list(...))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

setAs("SplitDataFrameList", "DataFrame", function(from) unlist(from))

setMethod("as.data.frame", "SplitDataFrameList",
          function(x, row.names=NULL, optional=FALSE, ...)
          {
            if (!(is.null(row.names) || is.character(row.names)))
              stop("'row.names'  must be NULL or a character vector")
            if (!missing(optional) || length(list(...)))
              warning("'optional' and arguments in '...' ignored")
            as.data.frame(as(x, "DataFrame"), row.names = row.names)
          })