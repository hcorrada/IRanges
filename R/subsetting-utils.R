### =========================================================================
### Subsetting utility functions
### -------------------------------------------------------------------------


### Returns an integer vector with values >= 1 and <= N, where N = length(x)
### if 'byrow=FALSE' and N = nrow(x) if 'byrow=TRUE'.
normalizeSingleBracketSubscript <- function(i, x, byrow=FALSE, exact=TRUE,
                                            allow.append=FALSE)
{
    if (!isTRUEorFALSE(byrow))
        stop("'byrow' must be TRUE or FALSE")
    if (!isTRUEorFALSE(exact))
        stop("'exact' must be TRUE or FALSE")
    if (!isTRUEorFALSE(allow.append))
        stop("'allow.append' must be TRUE or FALSE")
    if (byrow) {
        N <- nrow(x)
    } else {
        N <- length(x)
    }
    if (missing(i))
        return(seq_len(N))
    if (is.null(i))
        return(integer(0))
    if (is(i, "Rle")) {
        i <- as.vector(i)
    } else if (allow.append && is(i, "Ranges")) {
        i <- as.integer(i)
    }
    if (!is.atomic(i))
        stop("invalid subscript type")
    if (is.numeric(i)) {
        if (!is.integer(i))
            i <- as.integer(i)
        if (allow.append) {
            if (any(is.na(i)))
                stop("subscript contains NAs")
        } else {
            if (anyMissingOrOutside(i, upper=N))
                stop("subscript contains NAs or out of bounds indices")
        }
        nonzero_idx <- which(i != 0L)
        i <- i[nonzero_idx]
        if (length(i) == 0L)
            return(i)
        any_pos <- any(i > 0L)
        any_neg <- any(i < 0L)
        if (any_neg && any_pos)
            stop("cannot mix negative with positive indices")
        ## From here, indices are guaranteed to be either all positive or
        ## all negative.
        if (any_neg)
            return(seq_len(N)[i])
        return(i)
    }
    if (is.logical(i)) {
        if (anyMissing(i))
            stop("subscript contains NAs")
        li <- length(i)
        if (!allow.append && li > N) {
            if (any(i[(N+1L):li]))
                stop("subscript is a logical vector with out of bounds ",
                     "TRUE values")
            i <- i[seq_len(N)]
        }
        if (li < N)
            i <- rep(i, length.out=N)
        return(which(i))
    }
    if (is.character(i) || is.factor(i)) {
        if (byrow) {
            x_names <- rownames(x)
            what <- "rownames"
        } else {
            x_names <- names(x)
            what <- "names"
        }
        if (is.null(x_names)) {
            if (!allow.append)
                stop("cannot subset by character when ", what, " are NULL")
            return(N + seq_along(i))
        }
        if (exact) {
            i <- match(i, x_names, incomparables=c(NA_character_, ""))
        } else {
            i <- pmatch(i, x_names, duplicates.ok=TRUE)
        }
        if (allow.append) {
            na_idx <- which(is.na(i))
            i[na_idx] <- N + seq_along(na_idx)
            return(i)
        }
        if (anyMissing(i))
            stop("subscript contains invalid ", what)
        return(i)
    }
    stop("invalid subscript type")
}

### Supported types for 'i': single NA, numeric and character vectors only.
### Always returns a single integer. When called with 'error.if.nomatch=FALSE',
### returns an NA_integer_ if no match is found. Otherwise (the default),
### raises an error if no match is found so the returned integer is guaranteed
### to be a non-NA positive integer referring to a valid position in 'x'.
normalizeDoubleBracketSubscript <- function(i, x, exact=TRUE,
                                            error.if.nomatch=TRUE)
{
    if (!isTRUEorFALSE(exact))
        stop("'exact' must be TRUE or FALSE")
    if (!isTRUEorFALSE(error.if.nomatch))
        stop("'error.if.nomatch' must be TRUE or FALSE")
    if (missing(i))
        stop("subscript is missing")
    if (is.vector(i) && length(i) == 1L && is.na(i)) {
        if (error.if.nomatch)
            stop("subsetting by NA returns no match")
        return(NA_integer_)
    }
    if (!is.numeric(i) && !is.character(i))
        stop("invalid subscript type '", class(i), "'")
    if (length(i) < 1L)
        stop("attempt to extract less than one element")
    if (length(i) > 1L)
        stop("attempt to extract more than one element")
    if (is.numeric(i)) {
        if (!is.integer(i))
            i <- as.integer(i)
        if (i < 1L || length(x) < i)
            stop("subscript out of bounds")
        return(i)
    }
    ## 'i' is a character string
    x_names <- names(x)
    if (is.null(x_names)) {
        if (error.if.nomatch)
            stop("attempt to extract by name when elements have no names")
        return(NA_integer_)
    }
    #if (i == "")
    #    stop("invalid subscript \"\"")
    if (exact) {
        ans <- match(i, x_names, incomparables=c(NA_character_, ""))
    } else {
        ## Because 'i' has length 1, it doesn't matter whether we use
        ## 'duplicates.ok=FALSE' (the default) or 'duplicates.ok=TRUE' but
        ## the latter seems to be just a little bit faster.
        ans <- pmatch(i, x_names, duplicates.ok=TRUE)
    }
    if (is.na(ans) && error.if.nomatch)
        stop("subscript \"", i, "\" matches no name")
    ans
}

### Dispatch on the 2nd argument!
setGeneric("normalizeSingleBracketReplacementValue", signature="x",
    function(value, x, i)
        standardGeneric("normalizeSingleBracketReplacementValue")
)

### Default method.
setMethod("normalizeSingleBracketReplacementValue", "ANY",
    function(value, x)
    {
        if (is(value, class(x)))
            return(value)
        lv <- length(value)
        value <- try(as(value, class(x)), silent=TRUE)
        if (inherits(value, "try-error"))
            stop("'value' must be a ", class(x), " object (or coercible ",
                 "to a ", class(x), " object)")
        if (length(value) != lv)
            stop("coercing replacement value to ", class(x), "\n",
                 "  changed its length!\n",
                 "  Please do the explicit coercion ",
                 "yourself with something like:\n",
                 "    x[...] <- as(value, \"", class(x), "\")\n",
                 "  but first make sure this coercion does what you want.")
        value
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 2 internal generics to ease implementation of [ and [<- subsetting for
### new Vector subclasses.
###
### Most new Vector subclasses should only need to implement an "extractROWS"
### and a "replaceROWS" method to have "[" and "[<-" work out-of-the-box,
### respectively.
### Must support the following 'i' types: missing, Ranges and anything that
### can be handled by normalizeSingleBracketSubscript().
### For replaceROWS(), it's OK to assume that 'value' is "compatible" i.e.
### that it has gone thru normalizeSingleBracketReplacementValue().
### See "extractROWS" and "replaceROWS" methods for IRanges objects for an
### example.
###

setGeneric("extractROWS", signature="x",
    function(x, i) standardGeneric("extractROWS")
)

setGeneric("replaceROWS", signature="x",
    function(x, i, value) standardGeneric("replaceROWS")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 2 internal generics to ease implementation of [[ and [[<- subsetting for
### new List subclasses.
###

setGeneric("getListElement", signature="x",
    function(x, i, exact=TRUE) standardGeneric("getListElement")
)

setGeneric("setListElement", signature="x",
    function(x, i, value) standardGeneric("setListElement")
)

