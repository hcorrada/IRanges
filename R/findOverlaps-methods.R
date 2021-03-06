### =========================================================================
### findOverlaps (and related) methods
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### findOverlaps()
###
### Find objects in the query that overlap those in the subject.
###

setGeneric("findOverlaps", signature = c("query", "subject"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within", "equal"),
             select = c("all", "first", "last", "arbitrary"), ...)
        standardGeneric("findOverlaps")
)

###
### The next two functions are pre/post processing functions shared
### by the IntervalTree and IntervalForest versions of findOverlaps
###
### pre-process query (sort and adjust) and partition (if not NULL)
.preProcess_findOverlaps_query <- function(query, maxgap, minoverlap, partitioning=NULL) {
  res <- list()  
  query <- as(query, "IRanges")
  query_ord <- NULL
  res$origQuery <- query
  adjust <- (maxgap - minoverlap + 1L) * (width(query) > 0L)
  if (sum(adjust) > 0L)
    query <-
      resize(query, width(query) + 2L * adjust, fix = "center")
  res$unsortedQuery <- query

  if (is.null(partitioning)) {
    if (isNotSorted(start(query))) { ## query must be sorted
      query_ord <- sort.list(start(query), method = "quick",
                             na.last = NA)
      query <- query[query_ord]
    } else {
      query_ord <- seq_len(length(query))
    }
  } else {
    if(!is(partitioning, "PartitioningByEnd"))
      stop("invalid partitioning")

    # query be sorted within partition spaces 
    query_ord <- seq_len(length(query))

    isSorted <- TRUE
    for (i in seq(len=length(partitioning))) {
      curStart <- start(partitioning)[i]
      curEnd <- end(partitioning)[i]
      cur <- IRanges(curStart, curEnd)

      qStarts <- extractROWS(start(query), cur)
      if (isNotSorted(qStarts)) {
        isSorted <- FALSE
        ind <- extractROWS(query_ord, cur)
        query_ord <- replaceROWS(query_ord, cur,
                         ind[sort.list(qStarts, method="quick", na.last=NA)])
      }
    }
    if (!isSorted) {
      query <- query[query_ord]
    }
  }
  res$query <- query
  res$query_ord <- query_ord
  res
}

### post-process result based on overlap type and select
.postProcess_findOverlaps_result <- function(result, unsortedQuery, origQuery, subject, type, minoverlap, maxgap, origSelect)
{
  if (type != "any" || minoverlap > 1L) {
    if (!is(subject,"Ranges"))
        subject <- as(subject, "IRanges")

    m <- as.matrix(result)
    if (minoverlap > 1L) {

      r <- ranges(result, unsortedQuery, subject)
      m <- m[width(r) >= minoverlap, , drop=FALSE]
      ## unname() required because in case 'm' has only 1 row
      ## 'm[ , 1L]' and 'm[ , 2L]' will return a named atomic vector
      result@queryHits <- unname(m[ , 1L])
      result@subjectHits <- unname(m[ , 2L])
    }
    query <- origQuery
    
    filterMatrix <- function(fun)
      m[abs(fun(query)[m[,1L]] - fun(subject)[m[,2L]]) <= maxgap, ,
        drop=FALSE]
    
    if (type == "within") {
      r <- ranges(result, query, subject)
      m <- m[width(query)[m[,1L]] - width(r) <= maxgap, , drop=FALSE]
    } else if (type == "start") {
      m <- filterMatrix(start)
    } else if (type == "end") {
      m <- filterMatrix(end)
    } else if (type == "equal") {
      m <- filterMatrix(start)
      m <- filterMatrix(end)
    }
    
    if (origSelect != "all") {
      m <- m[!duplicated(m[,1L]), , drop=FALSE]
      result <- rep.int(NA_integer_, length(query))
      ## unname() required because in case 'm' has only 1 row
      ## 'm[,2L]' will return a named atomic vector
      result[m[,1L]] <- unname(m[,2L])
    } else {
      ## unname() required because in case 'm' has only 1 row
      ## 'm[ , 1L]' and 'm[ , 2L]' will return a named atomic vector
      result@queryHits <- unname(m[ , 1L])
      result@subjectHits <- unname(m[ , 2L])
    }
  }
  result
}

### findOverlaps method for IntervalTree
setMethod("findOverlaps", c("Ranges", "IntervalTree"),
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end", "within", "equal"),
                   select = c("all", "first", "last", "arbitrary"))
          {
            # verify inputs
            if (!isSingleNumber(maxgap) || maxgap < 0L)
              stop("'maxgap' must be a single, non-negative integer")
            if (!isSingleNumber(minoverlap) || minoverlap < 1L)
              stop("'minoverlap' must be a single, positive integer")
            type <- match.arg(type)
            select <- match.arg(select)


            origSelect <- select
            if (type != "any" || minoverlap > 1L)
              select <- "all"

            # preprocess query
            preprocRes <- .preProcess_findOverlaps_query(query, maxgap, minoverlap)
            origQuery <- preprocRes$origQuery
            unsortedQuery <- preprocRes$origQuery
            query <- preprocRes$query
            query_ord <- preprocRes$query_ord

            validObject(query)

            # make initial findOverlaps call
            fun <- paste("overlap_", select, sep = "")
            result <- IntervalTreeCall(subject, fun, query, query_ord)

            # postprocess results
            .postProcess_findOverlaps_result(result, unsortedQuery, origQuery, subject, type, minoverlap, maxgap, origSelect)
          })

setMethod("findOverlaps", c("Ranges", "Ranges"),
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end", "within", "equal"),
                   select = c("all", "first", "last", "arbitrary"))
          {
            findOverlaps(query, IntervalTree(subject), maxgap = maxgap,
                         minoverlap = minoverlap, type = match.arg(type),
                         select = match.arg(select))
          })

setMethod("findOverlaps", c("Vector", "missing"),
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end", "within", "equal"),
                   select = c("all", "first", "last", "arbitrary"),
                   ignoreSelf = FALSE, ignoreRedundant = FALSE)
          {
            select <- match.arg(select)
            result <- findOverlaps(query, query,
                                   maxgap = maxgap, minoverlap = minoverlap,
                                   type = type, select = "all")
            processSelfMatching(result, select, ignoreSelf, ignoreRedundant)
          })

setMethod("findOverlaps", c("integer", "Ranges"),
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end", "within", "equal"),
                   select = c("all", "first", "last", "arbitrary"))
          {
            findOverlaps(IRanges(query, query), subject, maxgap = maxgap,
                         minoverlap = minoverlap, type = match.arg(type),
                         select = match.arg(select))
          })

setMethod("findOverlaps", c("Views", "Views"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"))
    {
        findOverlaps(ranges(query), ranges(subject),
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=type, select=select)
    }
)

setMethod("findOverlaps", c("Views", "Vector"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"))
    {
        findOverlaps(ranges(query), subject,
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=type, select=select)
    }
)

setMethod("findOverlaps", c("Vector", "Views"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"))
    {
        findOverlaps(query, ranges(subject),
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=type, select=select)
    }
)

### findOverlaps method for IntervalForest
setMethod("findOverlaps", c("RangesList", "IntervalForest"),
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end", "within", "equal"),
                   select = c("all", "first", "last", "arbitrary"), drop=FALSE)
          {
            # verify inputs
            if (!isSingleNumber(maxgap) || maxgap < 0L)
              stop("'maxgap' must be a single, non-negative integer")
            if (!isSingleNumber(minoverlap) || minoverlap < 1L)
              stop("'minoverlap' must be a single, positive integer")

            type <- match.arg(type)
            select <- match.arg(select)
            origSelect <- select
            if (type != "any" || minoverlap > 1L)
              select <- "all"

            if (!is(query, "CompressedIRangesList"))
              query <- as(query, "CompressedIRangesList")

            validObject(query)

            queryList <- query
            query <- queryList@unlistData
            partitioning <- queryList@partitioning

            # preprocess query
            preprocRes <- .preProcess_findOverlaps_query(query, maxgap, minoverlap, partitioning)
            origQuery <- preprocRes$origQuery
            unsortedQuery <- preprocRes$unsortedQuery
            query <- preprocRes$query
            query_ord <- preprocRes$query_ord

            validObject(query)

            # match query partition to subject partition            
            partitionIndices <- match(names(partitioning), names(subject))
            
            # make initial findOverlaps call
            fun <- paste("overlap_", select, sep = "")
            result <- .IntervalForestCall(subject, fun, query, partitionIndices, elementLengths(partitioning), query_ord)

            # postprocess findOverlaps result
            res <- .postProcess_findOverlaps_result(result, unsortedQuery, origQuery, subject, type, minoverlap, maxgap, origSelect)

            # turn it into compressed list            
            if (origSelect == "all")
              return(CompressedHitsList(res, queryList))

            if (!drop) {
              return(newCompressedList0("CompressedIntegerList", unlistData=res, partitioning=partitioning))
            }
            res
          })

# might consider making this the following:
# setMethod("findOverlaps", c("RangesList", "RangesList"),
#           function(query, subject, maxgap = 0L, minoverlap = 1L,
#                    type = c("any", "start", "end", "within", "equal"),
#                    select = c("all", "first", "last", "arbitrary"),
#                    drop = FALSE)
#           {
#             findOverlaps(query, as(query, "IntervalForest"),
#               maxgap = maxgap, minoverlap = minoverlap,
#               type = match.arg(type), select = match.arg(select), drop = drop)
#           }
#   )

setMethod("findOverlaps", c("RangesList", "RangesList"),
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end", "within", "equal"),
                   select = c("all", "first", "last", "arbitrary"),
                   drop = FALSE)
          {
            type <- match.arg(type)
            select <- match.arg(select)

            query <- as.list(query)
            subject <- as.list(subject)
            origSubject <- subject
            if (!is.null(names(subject)) && !is.null(names(query))) {
              subject <- subject[names(query)]
              names(subject) <- names(query) # get rid of NA's in names
            } else {
              subject <- subject[seq_along(query)]
            }
            ## NULL's are introduced where they do not match
            ## We replace those with empty IRanges
            subject[sapply(subject, is.null)] <- IRanges()

            ans <- lapply(seq_len(length(subject)), function(i) {
              findOverlaps(query[[i]], subject[[i]], maxgap = maxgap,
                           minoverlap = minoverlap, type = type,
                           select = select)
            })
            names(ans) <- names(subject)
            if (select == "all") {
              ans <- HitsList(ans, origSubject)
            } else if (drop) {
              off <- head(c(0L, cumsum(sapply(origSubject, length))), -1)
              names(off) <- names(origSubject)
              if (is.null(names(ans)))
                off <- off[seq_along(ans)]
              else
                off <- off[names(ans)]
              ans <-
                unlist(ans, use.names=FALSE) +
                  rep.int(unname(off), sapply(ans, length))
            } else {
              ans <- IntegerList(ans)
            }
            ans
          })


setMethod("findOverlaps", c("ViewsList", "ViewsList"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             drop=FALSE)
    {
        findOverlaps(ranges(query), ranges(subject),
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=type, select=select, drop=drop)
    }
)

setMethod("findOverlaps", c("ViewsList", "Vector"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             drop=FALSE)
    {
        findOverlaps(ranges(query), subject,
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=type, select=select, drop=drop)
    }
)

setMethod("findOverlaps", c("Vector", "ViewsList"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             drop=FALSE)
    {
        findOverlaps(query, ranges(subject),
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=type, select=select, drop=drop)
    }
)

setMethod("findOverlaps", c("RangedData", "RangedData"),
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end", "within", "equal"),
                   select = c("all", "first", "last", "arbitrary"),
                   drop = FALSE)
          {
            findOverlaps(ranges(query), ranges(subject), maxgap = maxgap,
                         minoverlap = minoverlap, type = match.arg(type),
                         select = match.arg(select), drop = drop)
          })

setMethod("findOverlaps", c("RangedData", "RangesList"),
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end", "within", "equal"),
                   select = c("all", "first", "last", "arbitrary"),
                   drop = FALSE)
          {
            findOverlaps(ranges(query), subject, maxgap = maxgap,
                         minoverlap = minoverlap, type = match.arg(type),
                         select = match.arg(select), drop = drop)
          })

setMethod("findOverlaps", c("RangesList", "RangedData"),
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end", "within", "equal"),
                   select = c("all", "first", "last", "arbitrary"),
                   drop = FALSE)
          {
            findOverlaps(query, ranges(subject), maxgap = maxgap,
                         minoverlap = minoverlap, type = match.arg(type),
                         select = match.arg(select), drop = drop)
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### countOverlaps()
###

setGeneric("countOverlaps", signature = c("query", "subject"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within", "equal"), ...)
        standardGeneric("countOverlaps")
)

setMethod("countOverlaps", c("ANY", "Vector"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within", "equal"), ...)
    {
        counts <- queryHits(findOverlaps(query, subject, maxgap = maxgap,
                                         minoverlap = minoverlap, type = type, ...))
        structure(tabulate(counts, length(query)), names=names(query))

    }
)

setMethod("countOverlaps", c("ANY", "missing"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within", "equal"))
    {
        countOverlaps(query, query, maxgap = maxgap,
                      minoverlap = minoverlap, type = type)
    }
)

setMethod("countOverlaps", c("RangesList", "IntervalForest"),
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end", "within", "equal"), drop = FALSE, ...)
          {
            if (!is(query, "CompressedIRangesList"))
              query <- as("CompressedIRangesList", query)
            hits <- findOverlaps(query, subject, maxgap = maxgap,
                                 minoverlap = minoverlap, type = type, ...)
            res <- tabulate(queryHits(hits), queryLength(hits))
            new2("CompressedIntegerList", unlistData=res, partitioning = query@partitioning)
          })
  
setMethod("countOverlaps", c("RangesList", "RangesList"),
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end", "within", "equal"))
          {
              IntegerList(mapply(countOverlaps, query, subject,
                              MoreArgs = list(maxgap = maxgap,
                                      minoverlap = minoverlap,
                                      type = match.arg(type)),
                              SIMPLIFY = FALSE))
          })

setMethod("countOverlaps", c("ViewsList", "ViewsList"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"))
    {
         countOverlaps(ranges(query), ranges(subject),
                       maxgap=maxgap, minoverlap=minoverlap,
                       type=type)
    }
)

setMethod("countOverlaps", c("ViewsList", "Vector"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"))
    {
         countOverlaps(ranges(query), subject,
                       maxgap=maxgap, minoverlap=minoverlap,
                       type=type)
    }
)

setMethod("countOverlaps", c("Vector", "ViewsList"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"))
    {
         countOverlaps(query, ranges(subject),
                       maxgap=maxgap, minoverlap=minoverlap,
                       type=type)
    }
)

setMethod("countOverlaps", c("RangedData", "RangedData"),
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end", "within", "equal"))
          {
              countOverlaps(ranges(query), ranges(subject), maxgap = maxgap,
                            minoverlap = minoverlap, type = match.arg(type))
          })
setMethod("countOverlaps", c("RangedData", "RangesList"),
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end", "within", "equal"))
          {
              countOverlaps(ranges(query), subject, maxgap = maxgap,
                            minoverlap = minoverlap, type = match.arg(type))
          })
setMethod("countOverlaps", c("RangesList", "RangedData"),
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end", "within", "equal"))
          {
              countOverlaps(query, ranges(subject), maxgap = maxgap,
                            minoverlap = minoverlap, type = match.arg(type))
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### match() is defunct
###

setMethod("match", c("Views", "Views"),
    function(x, table, nomatch=NA_integer_, incomparables=NULL)
    {
        if (!identical(nomatch, NA_integer_))
            stop("'nomatch' arg is not supported")
        msg <- c("match() on Views objects is defunct.\n",
                 "Please use '",
                 "findOverlaps(ranges(x), ranges(table), select=\"first\")",
                 "' instead.")
        .Defunct(msg=msg)
    }
)

setMethod("match", c("Views", "Vector"),
    function(x, table, nomatch=NA_integer_, incomparables=NULL)
    {
        if (!identical(nomatch, NA_integer_))
            stop("'nomatch' arg is not supported")
        msg <- c("match() on a Views query is defunct.\n",
                 "Please use '",
                 "findOverlaps(ranges(x), table, select=\"first\")",
                 "' instead.")
        .Defunct(msg=msg)
    }
)

setMethod("match", c("Vector", "Views"),
    function(x, table, nomatch=NA_integer_, incomparables=NULL)
    {
        if (!identical(nomatch, NA_integer_))
            stop("'nomatch' arg is not supported")
        msg <- c("match() on a Views subject is defunct.\n",
                 "Please use '",
                 "findOverlaps(x, ranges(table), select=\"first\")",
                 "' instead.")
        .Defunct(msg=msg)
    }
)

setMethod("match", c("RangesList", "RangesList"),
          function(x, table, nomatch = NA_integer_, incomparables = NULL)
          {
            if (!identical(nomatch, NA_integer_))
                stop("'nomatch' arg is not supported")
            msg <- c("match() on RangesList objects is defunct.\n",
                     "Please use '",
                     "findOverlaps(x, table, select=\"first\", drop=TRUE)",
                     "' instead.")
            .Defunct(msg=msg)
          })

setMethod("match", c("ViewsList", "ViewsList"),
    function(x, table, nomatch=NA_integer_, incomparables=NULL)
    {
        if (!identical(nomatch, NA_integer_))
            stop("'nomatch' arg is not supported")
        msg <- c("match() on ViewsList objects is defunct.\n",
                 "Please use '",
                 "findOverlaps(ranges(x), ranges(table), select=\"first\", ",
                 "drop=TRUE)",
                 "' instead.")
        .Defunct(msg=msg)
    }
)

setMethod("match", c("ViewsList", "Vector"),
    function(x, table, nomatch=NA_integer_, incomparables=NULL)
    {
        if (!identical(nomatch, NA_integer_))
            stop("'nomatch' arg is not supported")
        msg <- c("match() on a ViewsList query is defunct.\n",
                 "Please use '",
                 "findOverlaps(ranges(x), table, select=\"first\", ",
                 "drop=TRUE)",
                 "' instead.")
        .Defunct(msg=msg)
    }
)

setMethod("match", c("Vector", "ViewsList"),
    function(x, table, nomatch=NA_integer_, incomparables=NULL)
    {
        if (!identical(nomatch, NA_integer_))
            stop("'nomatch' arg is not supported")
        msg <- c("match() on a ViewsList subject is defunct.\n",
                 "Please use '",
                 "findOverlaps(x, ranges(table), select=\"first\", ",
                 "drop=TRUE)",
                 "' instead.")
        .Defunct(msg=msg)
    }
)

setMethod("match", c("RangedData", "RangedData"),
    function(x, table, nomatch = NA_integer_, incomparables = NULL)
    {
        if (!identical(nomatch, NA_integer_))
            stop("'nomatch' arg is not supported")
        msg <- c("match() on RangedData objects is defunct.\n",
                 "Please use '",
                 "findOverlaps(ranges(x), ranges(table), select=\"first\", ",
                 "drop=TRUE)",
                 "' instead.")
        .Defunct(msg=msg)
    }
)

setMethod("match", c("RangedData", "RangesList"),
    function(x, table, nomatch = NA_integer_, incomparables = NULL)
    {
        if (!identical(nomatch, NA_integer_))
            stop("'nomatch' arg is not supported")
        msg <- c("match() on a RangedData query is defunct.\n",
                 "Please use '",
                 "findOverlaps(ranges(x), table, select=\"first\", ",
                 "drop=TRUE)",
                 "' instead.")
        .Defunct(msg=msg)
    }
)

setMethod("match", c("RangesList", "RangedData"),
    function(x, table, nomatch = NA_integer_, incomparables = NULL)
    {
        if (!identical(nomatch, NA_integer_))
            stop("'nomatch' arg is not supported")
        msg <- c("match() on a RangedData subject is defunct.\n",
                 "Please use '",
                 "findOverlaps(x, ranges(table), select=\"first\", ",
                 "drop=TRUE)",
                 "' instead.")
        .Defunct(msg=msg)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### overlapsAny()
###
### %in% is defunct. Replacing it with %over%.
###

### Same args and signature as countOverlaps() and subsetByOverlaps().
setGeneric("overlapsAny", signature=c("query", "subject"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"), ...)
        standardGeneric("overlapsAny")
)

setMethod("overlapsAny", c("Ranges", "Ranges"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"), ...)
    {
        !is.na(findOverlaps(query, subject, maxgap=maxgap,
                            minoverlap=minoverlap, type=type,
                            select="arbitrary", ...))
    }
)

setMethod("overlapsAny", c("Views", "Views"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"), ...)
    {
        overlapsAny(ranges(query), ranges(subject), maxgap=maxgap,
                    minoverlap=minoverlap, type=type, ...)
    }
)

setMethod("overlapsAny", c("Views", "Vector"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"), ...)
    {
        overlapsAny(ranges(query), subject, maxgap=maxgap,
                    minoverlap=minoverlap, type=type, ...)
    }
)

setMethod("overlapsAny", c("Vector", "Views"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"), ...)
    {
        overlapsAny(query, ranges(subject), maxgap=maxgap,
                    minoverlap=minoverlap, type=type, ...)
    }
)

setMethod("overlapsAny", c("RangesList", "RangesList"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"), ...)
    {
        query <- as.list(query)
        subject <- as.list(subject)
        if (!is.null(names(query)) && !is.null(names(subject))) {
            subject <- subject[names(query)]
            names(subject) <- names(query) # get rid of NA's in names
        } else {
            subject <- subject[seq_along(query)]
        }
        ## NULL's are introduced where they do not match
        ## We replace those with empty IRanges
        subject[sapply(subject, is.null)] <- IRanges()
        LogicalList(lapply(structure(seq_len(length(query)),
                                     names = names(query)),
                           function(i)
                               overlapsAny(query[[i]], subject[[i]],
                                           maxgap=maxgap,
                                           minoverlap=minoverlap,
                                           type=type, ...)))
    }
)

setMethod("overlapsAny", c("RangesList", "IntervalForest"),
      function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"), ...)
    {

      if (!is(query, "CompressedIRangesList"))
        query <- as(query, "CompressedIRangesList")

      !is.na(findOverlaps(query, subject, maxgap=maxgap,
                            minoverlap=minoverlap, type=type,
                            select="arbitrary", ...))
    }
)

setMethod("overlapsAny", c("ViewsList", "ViewsList"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"), ...)
    {
        overlapsAny(ranges(query), ranges(subject), maxgap=maxgap,
                    minoverlap=minoverlap, type=type, ...)
    }
)

setMethod("overlapsAny", c("ViewsList", "Vector"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"), ...)
    {
        overlapsAny(ranges(query), subject, maxgap=maxgap,
                    minoverlap=minoverlap, type=type, ...)
    }
)

setMethod("overlapsAny", c("Vector", "ViewsList"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"), ...)
    {
        overlapsAny(query, ranges(subject), maxgap=maxgap,
                    minoverlap=minoverlap, type=type, ...)
    }
)

setMethod("overlapsAny", c("RangedData", "RangedData"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"), ...)
    {
        overlapsAny(ranges(query), ranges(subject), maxgap=maxgap,
                    minoverlap=minoverlap, type=type, ...)
    }
)

setMethod("overlapsAny", c("RangedData", "RangesList"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"), ...)
    {
        overlapsAny(ranges(query), subject, maxgap=maxgap,
                    minoverlap=minoverlap, type=type, ...)
    }
)

setMethod("overlapsAny", c("RangesList", "RangedData"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"), ...)
    {
        overlapsAny(query, ranges(subject), maxgap=maxgap,
                    minoverlap=minoverlap, type=type, ...)
    }
)

### Convenience wrappers for the 2 most common use cases.
`%over%` <- function(query, subject) overlapsAny(query, subject)
`%within%` <- function(query, subject) overlapsAny(query, subject,
                                                   type="within")
`%outside%` <- function(query, subject) !overlapsAny(query, subject)

`.%in%.definition` <- function(x, table)
{
    msg <- c("%in% between a ", class(x), " and a ", class(table), " object ",
             "is defunct.\nPlease use 'query %over% subject' instead.")
    .Defunct(msg=msg)
}

.signatures <- list(
    c("Views", "Views"),
    c("Views", "Vector"),
    c("Vector", "Views"),
    c("RangesList", "RangesList"),
    c("ViewsList", "ViewsList"),
    c("ViewsList", "Vector"),
    c("Vector", "ViewsList"),
    c("RangedData", "RangedData"),
    c("RangedData", "RangesList"),
    c("RangesList", "RangedData")
)

setMethods("%in%", .signatures, `.%in%.definition`)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### subsetByOverlaps()
###

setGeneric("subsetByOverlaps", signature = c("query", "subject"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within", "equal"), ...)
        standardGeneric("subsetByOverlaps")
)

setMethod("subsetByOverlaps", c("Vector", "Vector"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within", "equal"), ...)
    {
        type <- match.arg(type)
        query[!is.na(findOverlaps(query, subject, maxgap = maxgap,
                                  minoverlap = minoverlap, type = type,
                                  select = "arbitrary", ...))]
    }
)

setMethod("subsetByOverlaps", c("RangedData", "RangedData"),
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end", "within", "equal"))
          {
              query[unlist(!is.na(findOverlaps(ranges(query), ranges(subject),
                                               maxgap = maxgap,
                                               minoverlap = minoverlap,
                                               type = match.arg(type),
                                               select = "arbitrary")),
                           use.names=FALSE),]
          })

setMethod("subsetByOverlaps", c("RangedData", "RangesList"),
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end", "within", "equal"))
          {
              query[unlist(!is.na(findOverlaps(ranges(query), subject,
                                               maxgap = maxgap,
                                               minoverlap = minoverlap,
                                               type = match.arg(type),
                                               select = "arbitrary")),
                           use.names=FALSE),]
          })

setMethod("subsetByOverlaps", c("RangesList", "RangedData"),
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                  type = c("any", "start", "end", "within", "equal"))
          {
              query[!is.na(findOverlaps(query, ranges(subject),
                                        maxgap = maxgap,
                                        minoverlap = minoverlap,
                                        type = match.arg(type),
                                        select = "arbitrary"))]
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "ranges" method for Hits objects
###
### Extracts the actual regions of intersection between the overlapping ranges.
### Not much value. Could be replaced by 1-liner:
###   pintersect(query[queryHits(x)], subject[subjectHits(x)])
###

setMethod("ranges", "Hits", function(x, query, subject) {
  if (!is(query, "Ranges") || length(query) != queryLength(x))
    stop("'query' must be a Ranges of length equal to number of queries")
  if (!is(subject, "Ranges") || length(subject) != subjectLength(x))
    stop("'subject' must be a Ranges of length equal to number of subjects")
  m <- as.matrix(x)
  qstart <- start(query)[m[,1L]]
  qend <- end(query)[m[,1L]]
  sstart <- start(subject)[m[,2L]]
  send <- end(subject)[m[,2L]]
  IRanges(pmax.int(qstart, sstart), pmin.int(send, qend))
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The only reason for defining the methods below is to prevent the default
### "findMatches" or "countMatches" methods to be called and return something
### wrong (and the reason they would return something wrong is because they
### are based on match() which does overlaps instead of equality).
### TODO: Remove these methods in BioC 2.14 when the "match" methods for all
### the signatures in '.signatures' are gone.

setMethods("findMatches", .signatures,
    function(x, table, select=c("all", "first", "last"), ...)
    {
        msg <- c("findMatches() between a ", class(x), " and a ",
                 class(table), " object is not supported")
        stop(msg)
    }
)

setMethods("countMatches", .signatures,
    function(x, table, ...)
    {
        msg <- c("countMatches() between a ", class(x), " and a ",
                 class(table), " object is not supported")
        stop(msg)
    }
)

