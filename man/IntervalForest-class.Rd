\name{IntervalForest-class}
\docType{class}
\alias{IntervalForest-class}

% constructor
\alias{IntervalForest}

% accessors
\alias{length,IntervalForest-method}
\alias{start,IntervalForest-method}
\alias{end,IntervalForest-method}

\title{Interval Search Forests}
\description{
  Efficiently perform overlap queries with a set of interval trees.
}
\details{
  A common type of query that arises when working with partitioned intervals is
  finding which intervals in one set overlap those in another. An
  efficient family of algorithms for answering such queries is known as
  the Interval Tree. The \code{IntervalForest} class stores a set of Interval Trees
  corresponding to intervals that are partitioned into disjoint sets. See the 
  \link{IntervalTree} class for the underlying Interval Tree data structure.

  A canonical example of partitioned ranges are \code{\link[GenomicRanges]{GenomicRanges}} 
  objects, where intervals are partitioned by their \code{seqnames}. See the 
  \link[GenomicRanges]{GIntervalTree} class to see the use of \code{IntervalForest} objects 
  in this case.

  The simplest approach for finding overlaps is to call the
  \code{\link{findOverlaps}} function on a \link{Ranges} object. 
  See the man page of \code{\link{findOverlaps-methods}}
  for how to use this and other related functions.
}

\section{Constructor}{
  \describe{
    \item{}{IntervalForest(ranges, partition): Creates an \code{IntervalForest} from the
      ranges in \code{ranges}, an object coercible to
      \code{IntervalTree}, such as an \code{\linkS4class{IRanges}} object, 
      and a \code{partition} factor such that an \code{IntervalTree} is constructed for 
      each level in the \code{partition} factor. \code{partition} can be vector or Rle object 
      (of class \code{factor}).
    }
  }
}


\section{Accessors}{
  \describe{
    \item{}{\code{length(x)}: Gets the number of ranges stored in the
      forest. This is a fast operation that does not bring the ranges into
      R.}
    \item{}{\code{start(x)}: Get the starts of the ranges.}
    \item{}{\code{end(x)}: Get the ends of the ranges.}
    \item{}{\code{x@partition}: The range partition factor as an \code{Rle} object.}
    \item{}{\code{levels(x)}: Get the levels of the range partition factor.}
  }
}

\author{Hector Corrada Bravo, Michael Lawrence}
\seealso{
  \code{\link{findOverlaps-methods}} for finding/counting interval overlaps between
  two partitioned "range-based" objects,
  \code{\linkS4class{Ranges}}, the parent of this class,
  \code{\linkS4class{Hits}}, set of hits between 2 vector-like objects,
  \code{\link[GenomicRanges]{GIntervalTree}}, which uses \code{IntervalForest} objects.
}
\examples{
  query <- IRanges(c(1, 4, 9), c(5, 7, 10))
  qpartition <- factor(c("a","a","b"))

  subject <- IRanges(c(2, 2, 10), c(2, 3, 12))
  spartition <- factor(c("a","a","b"))
  forest <- IntervalForest(subject, spartition)

  findOverlaps(query, forest, partition=qpartition)
}

\keyword{classes}
\keyword{methods}
