\name{findOverlaps-IntervalForest-methods}

\alias{findOverlaps-IntervalForest-methods}
\alias{findOverlaps,Ranges,IntervalForest-method}
\alias{countOverlaps,Ranges,IntervalForest-method}
\alias{overlapsAny,Ranges,IntervalForest-method}
\alias{subsetByOverlaps,Ranges,IntervalForest-method}

\title{Finding overlapping ranges in an Interval Forest}

\description{
  Various methods for finding/counting interval overlaps between two
  partitioned "range-based" objects: a query and a subject, where the subject is
  stored as an \link{IntervalForest} object.

  See \code{?`\link[GenomicRanges]{findOverlaps,GenomicRanges,GIntervalTree-method}`}
  in the GenomicRanges package for methods that operate on
  \link[GenomicRanges]{GRanges} or \link[GenomicRanges]{GRangesList}
  objects with subjects stored as \link[GenomicRanges]{GIntervalTree} objects which
  use \link{IntervalForest} objects as the backend.
}

\usage{
findOverlaps(query, subject, partition, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"), ...)

countOverlaps(query, subject, partition, maxgap=0L, minoverlap=1L,
              type=c("any", "start", "end", "within", "equal"), ...)

overlapsAny(query, subject, partition, maxgap=0L, minoverlap=1L,
            type=c("any", "start", "end", "within", "equal"), ...)

subsetByOverlaps(query, subject, partition, maxgap=0L, minoverlap=1L,
                 type=c("any", "start", "end", "within", "equal"), ...)
}

\arguments{
  \item{query}{
    A \link{Ranges} object. This method assumes the \code{query} ranges are partitioned as specified by the \code{partition} argument and overlaps are computed with ranges in \code{subject} in the same partition.

    If \code{query} is unsorted, it is sorted first (by partition and then by start position), so it is usually better to sort up-front, to avoid a sort with each \code{findOverlaps}
    call.
  }
  \item{subject}{
    An \link{IntervalForest} object.
  }

  \item{partition}{
    A \link{factor} stored either as a \code{vector} or \code{Rle} object. This argument partitions the \code{query} ranges such that overlaps are computed with ranges in \code{subject} with the same partition. This is converted to an \code{Rle} internally if necessary, so it is more efficient to pass this as an \code{Rle} object.
  }
  \item{maxgap, minoverlap}{
    Intervals with a separation of \code{maxgap} or less and a minimum
    of \code{minoverlap} overlapping positions, allowing for
    \code{maxgap}, are considered to be overlapping.  \code{maxgap}
    should be a scalar, non-negative, integer. \code{minoverlap}
    should be a scalar, positive integer.
  }
  \item{type}{
    By default, any overlap is accepted. By specifying the \code{type}
    parameter, one can select for specific types of overlap. The types
    correspond to operations in Allen's Interval Algebra (see
    references). If \code{type} is \code{start} or \code{end}, the
    intervals are required to have matching starts or ends,
    respectively. While this operation seems trivial, the naive
    implementation using \code{outer} would be much less
    efficient. Specifying \code{equal} as the type returns the
    intersection of the \code{start} and \code{end} matches. If
    \code{type} is \code{within}, the query interval must be wholly
    contained within the subject interval. Note that all matches must
    additionally satisfy the \code{minoverlap} constraint described above.

    The \code{maxgap} parameter has special meaning with the special
    overlap types. For \code{start}, \code{end}, and \code{equal},
    it specifies the maximum difference in the starts, ends or both,
    respectively. For \code{within}, it is the maximum amount by which
    the query may be wider than the subject.
  }
  \item{select}{
    When \code{select} is \code{"all"} (the default), the results are
    returned as a \link{Hits} object. When \code{select} is \code{"first"},
    \code{"last"}, or \code{"arbitrary"} the results are returned as an
    integer vector of length \code{query} containing the first, last,
    or arbitrary overlapping interval in \code{subject}, with \code{NA}
    indicating intervals that did not overlap any intervals in \code{subject}.
 
    If \code{select} is \code{"all"}, a \link{Hits} object is returned.
    For all other \code{select} the return value depends on the \code{drop}
    argument. When \code{select != "all" && !drop}, an \link{IntegerList}
    is returned, where each element of the result corresponds to a space in
    \code{query}. When\code{select != "all" && drop}, an integer vector is
    returned containing indices that are offset to align with the unlisted
    \code{query}.
  }
  \item{...}{
    Further arguments to be passed to or from other methods:
    \itemize{
      \item \code{drop}: All methods accept the \code{drop} argument
            (\code{FALSE} by default). See \code{select} argument above
            for the details.
      \item \code{ignoreSelf}, \code{ignoreRedundant}: When \code{subject}
            is omitted, the \code{ignoreSelf} and \code{ignoreRedundant}
            arguments (both \code{FALSE} by default) are allowed.
            See \code{query} and \code{subject} arguments above for the
            details.
    }
  }
  \item{x}{
    \link{Hits} object returned by \code{findOverlaps}.
  }
}

\details{
  A common type of query that arises when working with partitioned intervals is
  finding which intervals in one set overlap those in another. 

  The simplest approach for finding overlaps is to call the
  \code{\link{findOverlaps}} function on a partitioned \link{Ranges} object. A canonical example of partitioned ranges are \code{\link[GenomicRanges]{GenomicRanges}} object, where intervals are partitioned by their \code{seqnames}. See the \link[GenomicRanges]{GIntervalTree} to see the use \code{IntervalForest} objects in this case.

  The \code{IntervalForest} class stores a set of Interval Trees
  corresponding to intervals that are partitioned into disjoint sets. See the \link{IntervalTree} class for the underlying Interval Tree data structure which stores its ranges as a tree that is optimized for overlap queries.
  Thus, for repeated queries against the same subject, it is more
  efficient to create an \link{IntervalForest} once 
  using the \code{\link{IntervalForest}} constructor and then perform
  the queries against the \link{IntervalForest} instance.
}

\value{
  \code{findOverlaps} returns either a \link{Hits} object when
  \code{select="all"} (the default), or an integer vector when
  \code{select} is not \code{"all"}.

  \code{countOverlaps} returns the overlap hit count for each range
  in \code{query} using the specified \code{findOverlaps} parameters.

  \code{overlapsAny} finds the ranges in \code{query} that overlap any
  of the ranges in \code{subject}. It returns a logical vector of length equal to the number of
  ranges in \code{query}. 

  \code{subsetByOverlaps} returns the subset of \code{query} that
  has an overlap hit with a range in \code{subject} using the specified
  \code{findOverlaps} parameters.
}

\references{
  Allen's Interval Algebra:
  James F. Allen: Maintaining knowledge about temporal intervals. In:
  Communications of the ACM. 26/11/1983. ACM Press. S. 832-843, ISSN 0001-0782
}

\author{Hector Corrada Bravo, Michael Lawrence}

\seealso{
  \itemize{
    \item The \link{Hits} class for representing a set of hits between
          2 vector-like objects.

    \item \link[GenomicRanges]{findOverlaps,GenomicRanges,GIntervalTree-method}
          in the GenomicRanges package for methods that operate on
          \link[GenomicRanges]{GRanges} 
          objects using \code{IntervalForest} objects for querying.

    \item The \code{\link{IntervalTree}} class and constructor.

    \item The \link{Ranges} class.
  }
}

\examples{
query <- IRanges(c(1, 4, 9), c(5, 7, 10))
qpartition <- factor(c("a","a","b"))

subject <- IRanges(c(2, 2, 10), c(2, 3, 12))
spartition <- factor(c("a","a","b"))

forest <- IntervalForest(subject, spartition)

## ---------------------------------------------------------------------
## findOverlaps()
## ---------------------------------------------------------------------

## at most one hit per query
findOverlaps(query, forest, partition=qpartition, select = "first")
findOverlaps(query, forest, partition=qpartition, select = "last")
findOverlaps(query, forest, partition=qpartition, select = "arbitrary")

## overlap even if adjacent only
## (FIXME: the gap between 2 adjacent ranges should be still considered
## 0. So either we have an argument naming problem, or we should modify
## the handling of the 'maxgap' argument so that the user would need to
## specify maxgap = 0L to obtain the result below.)
findOverlaps(query, forest, partition=qpartition, maxgap = 1L)

## shortcut
findOverlaps(query, forest, partition=qpartition)

## alternative overlap types
query <- IRanges(c(1, 5, 3, 4), width=c(2, 2, 4, 6))
qpartition <- factor(c("a","a","b","b"))

subject <- IRanges(c(1, 3, 5, 6), width=c(4, 4, 5, 4))
spartition <- factor(c("a","a","b","b"))
forest <- IntervalForest(subject, spartition)

findOverlaps(query, forest, partition=qpartition, type = "start")
findOverlaps(query, forest, partition=qpartition, type = "start", maxgap = 1L)
findOverlaps(query, forest, partition=qpartition, type = "end", select = "first")
ov <- findOverlaps(query, forest, partition=qpartition, type = "within", maxgap = 1L)
ov

## ---------------------------------------------------------------------
## overlapsAny()
## ---------------------------------------------------------------------

overlapsAny(query, forest, partition=qpartition, type="start")
overlapsAny(query, forest, partition=qpartition, type="end")

## ---------------------------------------------------------------------
## subsetByOverlaps()
## ---------------------------------------------------------------------
subsetByOverlaps(query, forest, partition=qpartition)

## ---------------------------------------------------------------------
## countOverlaps()
## ---------------------------------------------------------------------
countOverlaps(query, forest, partition=qpartition)

## ---------------------------------------------------------------------
## "ranges" METHOD FOR Hits OBJECTS
## ---------------------------------------------------------------------

## extract the regions of intersection between the overlapping ranges
ranges(ov, query, forest)
}


\keyword{methods}
