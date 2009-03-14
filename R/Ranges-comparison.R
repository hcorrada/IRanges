### =========================================================================
### Comparison of Ranges objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Equality and related methods.
###

setMethod("==", signature(e1="Ranges", e2="Ranges"),
    function(e1, e2)
    {
        (start(e1) == start(e2)) & (width(e1) == width(e2))
    }
)

setMethod("!=", signature(e1="Ranges", e2="Ranges"),
    function(e1, e2)
    {
        ## Should be slightly more efficient than '!(e1 == e2)', at least in
        ## theory.
        (start(e1) != start(e2)) | (width(e1) != width(e2))
    }
)

### Need to explicitly define this generic otherwise the implicit generic in
### package "base" would dispatch on (x, incomparables).
setGeneric("duplicated", signature="x",
    function(x, incomparables=FALSE, ...) standardGeneric("duplicated")
)

### Note that this default method is very inefficient so efficient methods for
### the Ranges subclasses need to be implemented.
setMethod("duplicated", "Ranges",
    function(x, incomparables=FALSE, fromLast=FALSE, ...)
    {
        if (!identical(incomparables, FALSE))
            stop("\"duplicated\" method for Ranges objects only accepts 'incomparables=FALSE'")
        duplicated(data.frame(start=start(x),
                              width=width(x),
                              check.names=FALSE,
                              stringsAsFactors=FALSE),
                   fromLast=fromLast)
    }
)

### Need to explicitly define this generic otherwise the implicit generic in
### package "base" would dispatch on (x, incomparables).
setGeneric("unique", signature="x",
    function(x, incomparables=FALSE, ...) standardGeneric("unique")
)

### Relies on a "[" method for 'x'.
setMethod("unique", "Ranges",
    function(x, incomparables=FALSE, fromLast=FALSE, ...)
    {
        x[-which(duplicated(x, incomparables=incomparables, fromLast=fromLast))]
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Ordering and related methods.
###
### Ranges are ordered by starting position first and then by width.
### This way, the space of ranges is totally ordered.
### The "order", "sort" and "rank" methods for Ranges objects are consistent
### with this order.
###

setMethod("<=", signature(e1="Ranges", e2="Ranges"),
    function(e1, e2)
    {
         (start(e1) < start(e2)) | ((start(e1) == start(e2)) & (width(e1) <= width(e2)))
    }
)

setMethod(">=", signature(e1="Ranges", e2="Ranges"),
    function(e1, e2) {e2 <= e1}
)

setMethod("<", signature(e1="Ranges", e2="Ranges"),
    function(e1, e2)
    {
         (start(e1) < start(e2)) | ((start(e1) == start(e2)) & (width(e1) < width(e2)))
    }
)

setMethod(">", signature(e1="Ranges", e2="Ranges"),
    function(e1, e2) {e2 < e1}
)

### Need to explicitly define this generic otherwise the implicit generic in
### package "base" would dispatch on (na.last, decreasing).
### Note that dispatching on ... is supported starting with R 2.8.0 only.
setGeneric("order", signature="...",
    function(..., na.last=TRUE, decreasing=FALSE) standardGeneric("order")
)

setMethod("order", "Ranges",
    function(..., na.last=TRUE, decreasing=FALSE)
    {
        ## all arguments in '...' are guaranteed to be Ranges objects
        args <- list(...)
        order_args <- vector("list", 2L*length(args))
        order_args[2L*seq_len(length(args)) - 1L] <- lapply(args, start)
        order_args[2L*seq_len(length(args))] <- lapply(args, end)
        do.call(order, c(order_args, na.last=na.last, decreasing=decreasing))
    }
)

### Need to explicitly define this generic otherwise the implicit generic in
### package "base" would dispatch on (x, decreasing).
setGeneric("sort", signature="x",
    function(x, decreasing=FALSE, ...) standardGeneric("sort")
)

### Relies on a "[" method for 'x'.
setMethod("sort", "Ranges",
    function(x, decreasing=FALSE, ...) 
    {
        if (!isTRUEorFALSE(decreasing))
            stop("'decreasing' must be TRUE or FALSE")
        x[order(x, decreasing=decreasing)]
    }
)

### Need to explicitly define this generic otherwise the implicit generic in
### package "base" would dispatch on (x, na.last, ties.method).
setGeneric("rank", signature="x",
    function(x, na.last=TRUE, ties.method=c("average", "first", "random", "max", "min"))
        standardGeneric("rank")
)

setMethod("rank", "Ranges",
    function(x, na.last=TRUE, ties.method=c("average", "first", "random", "max", "min"))
    {
        if (!identical(ties.method, "first"))
            stop("only 'ties.method=\"first\"' is supported when ranking ranges")
        ox <- order(x)
        ## 'ans' is the reverse permutation of 'ox'
        ans <- integer(length(ox))
        ans[ox] <- seq_len(length(ox))
        ans
    }
)
