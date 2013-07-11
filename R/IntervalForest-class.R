### =========================================================================
### IntervalForest objects
### -------------------------------------------------------------------------

setClass("IntervalForest",
         representation(ptr = "externalptr", mode = "character", partition="Rle"),
         contains = "Ranges")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("length", "IntervalForest", function(x) .IntervalForestCall(x, "length"))
setMethod("start", "IntervalForest", function(x) .IntervalForestCall(x, "start"))
setMethod("end", "IntervalForest", function(x) .IntervalForestCall(x, "end"))
setMethod("levels", "IntervalForest", function(x) levels(x@partition))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

IntervalForest <- function(ranges, partition) {
  if (is(partition, "Rle")) {
    if (!is.factor(runValue(partition))) {
      stop("'partition' must be a 'factor' Rle or 'factor'")
    }
  } else {
    if (!is.factor(partition)) {
      stop("'partition' must be a 'factor' Rle or 'factor'")
    }
    partition <- Rle(partition)
  }
  
  validObject(ranges)
  validObject(partition)
  
  
  npartitions <- nlevels(partition)
  levels <- levels(partition)
  partitionIndices <- as.integer(partition)
  
  ptr <- .Call2("IntegerIntervalForest_new", ranges, partitionIndices, npartitions, PACKAGE="IRanges")
  new2("IntervalForest", ptr = ptr, mode="integer", partition=partition, check=FALSE)
}


### - - - - 
### Subsetting
###

setMethod("[", "IntervalForest",
          function(x, i, j, ..., drop=TRUE) {
            if (!missing(j) || length(list(...)) > 0L)
              stop("invalid subsetting")
            newRanges <- callGeneric(as(x, "IRanges"), i=i, ...)
            newPartition <- callGeneric(x@partition, i=i, ...)
            IntervalForest(newRanges, newPartition)
          })

### - - - -
### updating
### this allows support for shift, narrow, etc. methods
### - - - - 

setMethod("update", "IntervalForest",
    function(object, ...) 
        IntervalForest(update(as(object, "IRanges"), ...), object@partition))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Low-level utilities
###

.IntervalForestCall <- function(object, fun, ...) {
  validObject(object)
  fun <- paste("IntervalForest", fun, sep = "_")
  if (object@mode == "integer") {
    fun <- paste("Integer", fun, sep = "")
    .Call2(fun, object@ptr, ..., PACKAGE="IRanges")
  } else stop("unknown interval forest mode: ", object@mode)
}

## not for exporting, just a debugging utility
IntervalForestDump <- function(object) {
  cat("IntervalForest, levels: ", levels(object), "\n")
  .IntervalForestCall(object, "dump")
}