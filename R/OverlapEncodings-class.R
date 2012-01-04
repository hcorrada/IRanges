### =========================================================================
### OverlapEncodings objects
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compare()
###
### TODO: This should probably go somewhere else e.g. in Ranges-comparison.R
###

setGeneric("compare", function(x, y) standardGeneric("compare"))

### "Parallel" generalized comparison of 2 Ranges objects.
###   > x <- IRanges(1:11, width=4)
###   > y <- IRanges(6, 9)
###   > compare(x, y)
###    [1] -6 -5 -4 -4 -4  0  4  4  4  5  6
###   > compare(IRanges(4:6, width=6), y)
###   [1] -3 -2  1
###   > compare(IRanges(6:8, width=2), y)
###   [1] -1  2  3
###   > compare(x, y) < 0  # equivalent to x < y
###   > compare(x, y) == 0  # equivalent to x == y
###   > compare(x, y) > 0  # equivalent to x > y
### TODO: Seems like using compare() to implement "==", "!=", "<=", ">=",
### "<" and ">" methods for Ranges objects would make them slightly faster
### (between 1.5x and 2.5x) and also slightly more memory efficient.
setMethod("compare", c("Ranges", "Ranges"),
    function(x, y)
    {
        .Call2("Ranges_compare",
               start(x), width(x), start(y), width(y),
               PACKAGE="IRanges")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The OverlapEncodings class.
###
### The OverlapEncodings class is a container for storing a vector of OVM's.
### An OVM ("overlaps matrix") is an M x N matrix of 1-letter codes
### representing all the range-to-range overlaps between a query with M ranges
### and a subject with N ranges.
###

### Slots:
###   o Loffset ("left offset"): nb of cols on the left of the OVM that contain
###     only "m"'s.
###   o Roffset ("right offset"): nb of cols on the right of the OVM that
###     contain only "a"'s.
###   o encoding: linear sequence of symbols representing the trimmed OVM (i.e.
###     after removing Loffset cols on the left and Roffset cols on the right).
setClass("OverlapEncodings",
    contains="Vector",
    representation(
        Loffset="integer",    # no NAs, >= 0
        Roffset="integer",    # no NAs, >= 0
        encoding="character"  # no NAs
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### encodeOverlaps()
###

### Low-level utility.
###   > query <- IRanges(c(7, 15, 22), c(9, 19, 23))
###   > subject <- IRanges(c(1, 4, 15, 22), c(2, 9, 19, 25))
###   > encodeOverlaps1(query, subject, sparse.output=FALSE)
###   [1] "mjaa" "mmga" "mmmf"
###   > encodeOverlaps1(query, subject)  # Type III encoding
###   [1] ":jmm:agm:aaf:"
### TODO: Do we really need this? Same could be achieved with
### 'encodeOverlaps(IRangesList(query), IRangesList(subject))' except that
### with encodeOverlaps1() we can specify 'query.space' and 'subject.space'
### and use the 'sparse.output' and 'as.raw' flags to control the format of
### the output.
### Also the C code behind encodeOverlaps1() is at the heart of all the
### "encodeOverlaps" methods. Playing with encodeOverlaps1() with different
### inputs and combinations of the 'sparse.output' and 'as.raw' flags is
### educative and allows testing all parts of this C code.
encodeOverlaps1 <- function(query, subject,
                            query.space=NULL, subject.space=NULL,
                            sparse.output=TRUE, as.raw=FALSE)
{
    if (!isTRUEorFALSE(sparse.output))
        stop("'sparse.output' must be TRUE or FALSE")
    if (!isTRUEorFALSE(as.raw))
        stop("'as.raw' must be TRUE or FALSE")
    .Call2("encode_overlaps1",
           start(query), width(query), query.space,
           start(subject), width(subject), subject.space,
           sparse.output, as.raw,
           PACKAGE="IRanges")
}

### TODO: Put this in the (upcoming) man page for encodeOverlaps().
### A simple (but inefficient) implementation of the "findOverlaps" method for
### Ranges objects. Complexity and memory usage is M x N where M and N are the
### lengths of 'query' and 'subject', respectively.
findRangesOverlaps <- function(query, subject)
{
    ## WARNING: When using sparse.output=FALSE and as.raw=TRUE, the returned
    ## raw matrix is transposed!
    ocodes <- encodeOverlaps1(query, subject, sparse.output=FALSE, as.raw=TRUE)
    offsets <- which(charToRaw("c") <= ocodes & ocodes <= charToRaw("k")) - 1L
    q_hits <- offsets %/% nrow(ocodes) + 1L
    s_hits <- offsets %% nrow(ocodes) + 1L
    cbind(query=q_hits, subject=s_hits)
}

RangesList_encodeOverlaps <- function(query.starts, query.widths,
                                      subject.starts, subject.widths,
                                      query.spaces=NULL, subject.spaces=NULL)
{
    .Call2("RangesList_encode_overlaps",
           query.starts, query.widths, query.spaces,
           subject.starts, subject.widths, subject.spaces,
           PACKAGE="IRanges")
}

setGeneric("encodeOverlaps", signature=c("query", "subject"),
    function(query, subject) standardGeneric("encodeOverlaps")
)

### "Parallel" overlap encoding between 2 RangesList objects.
### Between many reads and one transcript.
###   > read1 <- IRanges(c(7, 15, 22), c(9, 19, 23))
###   > read2 <- IRanges(c(5, 15), c(9, 17))
###   > read3 <- IRanges(c(16, 22), c(19, 24))
###   > query <- IRangesList(read1, read2, read3)
###   > tx <- IRanges(c(1, 4, 15, 22), c(2, 9, 19, 25))
###   > subject <- IRangesList(tx)
###   > ocodes <- encodeOverlaps(query, subject)
###   > ocodes
###   [1] ":jmm:agm:aaf:" ":jm:af:"       ":jm:af:"
### Reads compatible with transcript 'tx':
###   ## Regex to use for reads with no gaps:
###   > pattern0 <- ":[fgij]:"
###   ## Regex to use for reads with 1 gap:
###   > pattern1 <- ":[jg].:.[gf]:"
###   ## Regex to use for reads with 2 gaps:
###   > pattern2 <- ":[jg]..:.g.:..[gf]:"
###   ## Regex to use for reads with up to 2 gaps:
###   > pattern012 <- ":([fgij]|[jg].:.[gf]|[jg]..:.g.:..[gf]):"
###   > grep(pattern012, ocodes)
###   [1] 1 2 3
### All the reads are compatible with this transcript!
setMethod("encodeOverlaps", c("RangesList", "RangesList"),
    function(query, subject)
    {
        RangesList_encodeOverlaps(as.list(start(query)),
                                  as.list(width(query)),
                                  as.list(start(subject)),
                                  as.list(width(subject)))
    }
)

setMethod("encodeOverlaps", c("RangesList", "Ranges"),
    function(query, subject)
    {
        RangesList_encodeOverlaps(as.list(start(query)),
                                  as.list(width(query)),
                                  as.list(start(subject)),
                                  as.list(width(subject)))
    }
)

setMethod("encodeOverlaps", c("Ranges", "RangesList"),
    function(query, subject)
    {
        RangesList_encodeOverlaps(as.list(start(query)),
                                  as.list(width(query)),
                                  as.list(start(subject)),
                                  as.list(width(subject)))
    }
)

###   > query <- IRanges(1:11, width=4)
###   > subject <- IRanges(6, 9)
###   > encodeOverlaps(query, subject)
###    [1] "a" "b" "c" "c" "c" "g" "k" "k" "k" "l" "m"
###   > encodeOverlaps(IRanges(4:6, width=6), subject)
###   [1] "d" "e" "h"
###   > encodeOverlaps(IRanges(6:8, width=2), subject)
###   [1] "f" "i" "j"
setMethod("encodeOverlaps", c("Ranges", "Ranges"),
    function(query, subject)
    {
        ### TODO: Add an extra arg to compare() to let the user choose the
        ### type of output i.e. numeric or 1-letter codes.
        ocodes <- compare(query, subject)
        safeExplode(rawToChar(as.raw(as.integer(charToRaw("g")) + ocodes)))
    }
)
