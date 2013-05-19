### =========================================================================
### IntervalForest objects
### -------------------------------------------------------------------------

setClass("IntervalForest",
         representation(ptr = "externalptr", mode = "character"),
         contains = "Ranges")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("length", "IntervalForest", function(x) .IntervalForestCall(x, "length"))

setMethod("start", "IntervalForest", function(x) .IntervalForestCall(x, "start"))
setMethod("end", "IntervalForest", function(x) .IntervalForestCall(x, "end"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

IntervalForest <- function(ranges, partition) {
  validObject(ranges)
  validObject(partition)
  
  if (!is(partition, "factor")) {
    stop("partition must be of class 'factor'")
  }
  npartitions <- nlevels(partition)
  partition <- as.integer(partition)
  
  ptr <- .Call2("IntegerIntervalForest_new", ranges, partition, npartitions, PACKAGE="IRanges")
  new2("IntervalForest", ptr = ptr, mode="integer", check=FALSE)
}


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
