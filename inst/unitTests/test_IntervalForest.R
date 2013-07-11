test_IntervalForest_construction <- function() {
  query <- IRanges(c(1, 3, 9), c(5, 7, 10))
  qpartition <- factor(c("a","a","b"))
  
  subject <- IRanges(c(2, 10), c(2, 12))
  spartition <- factor(c("a","b"))
  
  tree <- IntervalForest(subject, spartition)
  checkTrue(validObject(tree))
  checkIdentical(length(tree), 2L);
  checkIdentical(levels(tree), c("a","b"))
  
  tree <- IntervalForest(IRanges(), factor())
  checkTrue(validObject(tree))

  tree <- IntervalForest(IRanges(1, 0), factor(c("a")))
  checkIdentical(start(tree), 1L)
  checkIdentical(levels(tree), c("a"))
  checkTrue(validObject(tree))

  tree <- IntervalForest(IRanges(c(1, 1), c(1, 0)), factor(c("a","b")))
  checkIdentical(width(tree), c(1L, 0L))
  checkIdentical(levels(tree), c("a","b"))
  checkTrue(validObject(tree))

  tree <- IntervalForest(IRanges(1:10,width=1), factor(rep("a",10)))
  checkIdentical(start(tree), as.integer(1:10))
  checkIdentical(width(tree), rep(1L,10))
  checkTrue(validObject(tree))
  
  checkException(IntervalForest(), silent = TRUE)
  checkException(IntervalForest(subject, query, qpartition), silent = TRUE)
  checkException(IntervalForest(NULL, NULL), silent = TRUE)
}

test_IntervalForest_findOverlaps <- function() {
  ## a .....
  ## b    ....
  ## a        ..
  ## a  x
  ## b  xx
  ## a         xxx
  query <- IRanges(c(1, 4, 9), c(5, 7, 10))
  qpartition <- factor(c("a","b","a"))
  
  subject <- IRanges(c(2, 2, 10), c(2, 3, 12))
  spartition <- factor(c("a","a","b"))
  
  tree <- IntervalForest(subject, spartition)

  result <- findOverlaps(query, tree, partition=qpartition, select = "first")
  checkIdentical(result, c(1L, NA, NA))
  
  query <- IRanges(c(1, 4, 9), c(5, 7, 10))
  qpartition <- factor(c("a","a","b"))
  
  result <- findOverlaps(query, tree, partition=qpartition, select = "first")
  checkIdentical(result, c(1L, NA, 3L))
  
  result <- findOverlaps(query, tree, partition=qpartition, select = "last")
  checkIdentical(result, c(2L, NA, 3L))
  
  result <- findOverlaps(query, tree, partition=qpartition, select = "arbitrary")
  checkIdentical(result, c(2L, NA, 3L))

  checkOverlap <- function(a, q, s, r, c) {
    mat <- cbind(queryHits = as.integer(q), subjectHits = as.integer(s))
    checkIdentical(as.matrix(a), mat)
    checkIdentical(queryLength(a), as.integer(r))
    checkIdentical(subjectLength(a), as.integer(c))
  }

   result <- findOverlaps(query, tree, partition=qpartition)
   checkOverlap(result, c(1, 1, 3), c(1, 2, 3), 3, 3)
 
  ## with 'maxgap'
  result <- findOverlaps(query, tree, 1, partition=qpartition)
  checkOverlap(result, c(1, 1, 2, 3), c(1, 2, 2, 3), 3, 3)

  ## with 'minoverlap'
  result <- findOverlaps(query, tree, minoverlap = 3L, partition=qpartition)
  checkOverlap(result, integer(0), integer(0), 3, 3)
  result <- findOverlaps(query, tree, minoverlap = 2L, partition=qpartition)
  checkOverlap(result, 1, 2, 3, 3)
  result <- findOverlaps(query, tree, minoverlap = 2L, select = "first", partition=qpartition)
  checkIdentical(result, c(2L, NA, NA))
  result <- findOverlaps(query, tree, minoverlap = 2L, select = "last", partition=qpartition)
  checkIdentical(result, c(2L, NA, NA))
  result <- findOverlaps(query, tree, minoverlap = 2L, select = "arbitrary", partition=qpartition)
  checkIdentical(result, c(2L, NA, NA))

  ## empty query range
  #subject <- IRanges(c(2, 2, 10), c(2, 3, 12))
  #spartition <- factor(c("a","a","b"))
  query <- IRanges(c(1, 4, 9, 10), c(5, 7, 10, 9))
  qpartition <- factor(c("a","a","b","b"))
  result <- findOverlaps(query, tree, partition=qpartition)
  checkOverlap(result, c(1, 1, 3), c(1, 2, 3), 4, 3)

  ## empty subject range
  subject <- IRanges(c(2, 2, 2, 10), c(2, 1, 3, 12))
  spartition <- factor(c("a","a","a","b"))
  tree <- IntervalForest(subject, spartition)
  result <- findOverlaps(query, tree, partition=qpartition)
  checkOverlap(result, c(1, 1, 3), c(1, 3, 4), 4, 4)

  ## .....
  ##    ....
  ##         ..
  ##  xxxx
  ##  xxx
  query <- IRanges(c(1, 4, 9), c(5, 7, 10))
  qpartition <- factor(c("a","b","a"))
  
  subject <- IRanges(c(2, 2), c(5, 4))
  spartition <- factor(c("a","b"))
  
  tree <- IntervalForest(subject, spartition)
  result <- findOverlaps(query, tree, partition=qpartition)
  checkOverlap(result, c(1, 2), c(1, 2), 3, 2)

  
  ## check case of identical subjects
  ## .....
  ##    .....
  ##         ..
  ##  xxxx
  ##  xxxx
  ##      xx
  ##      xxx
  ##      xx
  query <- IRanges(c(1, 4, 9), c(5, 7, 10))
  qpartition <- factor(c("a","b","a"))
  
  subject <- IRanges(c(2, 2, 6, 6, 6), c(5, 5, 7, 8, 7))  
  spartition <- factor(c("a","a","b","b","b"))
  
  tree <- IntervalForest(subject, spartition)
  result <- findOverlaps(query, tree, partition=qpartition)
  checkOverlap(result, c(1, 1, 2, 2, 2), c(1, 2, 3, 4, 5), 3, 5)

  # on unsorted query
  query <- IRanges(c(10, 5, 3, 7, 9), c(15, 7, 7, 10, 12))
  qpartition <- factor(c("a","b","b","a","a"))
  result <- findOverlaps(query, tree, partition=qpartition)
  checkOverlap(result, c(2, 2, 2, 3, 3, 3), c(3, 4, 5, 3, 4, 5), 5, 5)

  # query with partition level not in subject
  query <- IRanges(c(10, 5, 3, 7, 9), c(15, 7, 7, 10, 12))
  qpartition <- factor(c("a","b","c","a","a"))

  subject <- IRanges(c(2, 2, 6, 6, 6), c(5, 5, 7, 8, 7))  
  spartition <- factor(c("b","b","b","b","b"))

  tree <- IntervalForest(subject, spartition)
  result <- findOverlaps(query, tree, partition=qpartition)
  checkOverlap(result, c(2, 2, 2, 2, 2), c(1, 2, 3, 4, 5), 5, 5)

#   subject <- IRanges(c(1, 6, 13), c(4, 9, 14)) # single points (this doesn't work since findOverlaps-integer,Ranges doesn't pass ...)
#   spartition <- factor(c("a","b","c"))
#   subject <- IntervalForest(subject, spartition)
#   
#   checkIdentical(findOverlaps(c(3L, 7L, 10L), subject, partition=factor(c("a","a","b")), select = "first"),
#                  c(1L, 2L, NA))
#   checkIdentical(findOverlaps(c(3L, 7L, 10L), subject, partition=factor(c("a","a","b")), select = "last"),
#                  c(1L, 2L, NA))
#   checkIdentical(findOverlaps(c(3L, 7L, 10L), subject, partition=factor(c("a","a","b")), select = "arbitrary"),
#                  c(1L, 2L, NA))
#   checkIdentical(findOverlaps(IRanges(c(2,1),c(3,4)), subject, partition=factor(c("a","a"))),
#                  new("Hits",
#                      queryHits = 1:2, subjectHits = c(1L,1L),
#                      queryLength = 2L, subjectLength = 3L))

  ## check other types of matching

  ## ..
  ##     ..
  ##   ....  
  ##    ......
  ## xxxx
  ##   xxxx
  ##     xxxxx
  ##      xxxx

  query <- IRanges(c(1, 5, 3, 4), width=c(2, 2, 4, 6))
  qpartition <- factor(c("a","a","a","a"))
  
  subject <- IRanges(c(1, 3, 5, 6), width=c(4, 4, 5, 4))
  spartition <- factor(c("a","a","a","a"))
  tree <- IntervalForest(subject, spartition)

  ## 'start'
  result <- findOverlaps(query, tree, type = "start", partition=qpartition)
  checkOverlap(result, c(1, 2, 3), c(1, 3, 2), 4, 4)

  ## non-zero maxgap
  result <- findOverlaps(query, tree, type = "start", maxgap = 1L, partition=qpartition)
  checkOverlap(result, c(1, 2, 2, 3, 4, 4), c(1, 3, 4, 2, 2, 3), 4, 4)

  ## minoverlap > 1L
  result <- findOverlaps(query, tree, type = "start", minoverlap = 3L, partition=qpartition)
  checkOverlap(result, 3, 2, 4, 4)
  
  ## combine minoverlap and maxgap
  result <- findOverlaps(query, tree, type = "start", maxgap = 1L,
                         minoverlap = 3L, partition=qpartition)
  checkOverlap(result, c(3, 4, 4), c(2, 2, 3), 4, 4)
  
  ## 'end'
  result <- findOverlaps(query, tree, type = "end", partition=qpartition)
  checkOverlap(result, c(2, 3, 4, 4), c(2, 2, 3, 4), 4, 4)

#   ## ensure inverse is same as transpose
#   inverse <- findOverlaps(subject, query, type = "end")
#   tr <- as.matrix(t(result))
#   checkIdentical(as.matrix(inverse), tr[order(tr[,1]),])

  ## select = "first"
  result <- findOverlaps(query, tree, type = "end", select = "first",partition=qpartition)
  checkIdentical(result, c(NA, 2L, 2L, 3L))  

  ## 'within'
  result <- findOverlaps(query, tree, type = "within", partition=qpartition)
  checkOverlap(result, c(1, 2, 2, 3), c(1, 2, 3, 2), 4, 4)  

  result <- findOverlaps(query, tree, type = "within", maxgap = 1L, partition=qpartition)
  checkOverlap(result, c(1, 2, 2, 2, 3, 4), c(1, 2, 3, 4, 2, 3), 4, 4)  
  
  ## 'equal'
  result <- findOverlaps(query, tree, type = "equal", partition=qpartition)
  checkOverlap(result, 3, 2, 4, 4)  
# 
#   checkException(findOverlaps(query, NULL), silent = TRUE)
#   checkException(findOverlaps(NULL, query), silent = TRUE)
}

test_IntervalForest_asRanges <- function() {
  ranges <- IRanges(c(1, 4, 9), c(5, 7, 10))
  partition <- factor(c("a","b","a"))
  
  tree <- IntervalForest(ranges, partition)
  checkIdentical(as(tree, "IRanges"), ranges)
  
  ranges <- IRanges()
  partition <- factor()
  tree <- IntervalForest(ranges, partition)
  checkIdentical(as(tree, "IRanges"), ranges)
}

test_IntervalForest_subset <- function() {
  ranges <- IRanges(c(1, 4, 9), c(5, 7, 10))
  partition <- Rle(factor(c("a","b","a")))
  
  tree <- IntervalForest(ranges, partition)
  checkIdentical(as(tree, "IRanges"), ranges)
  
  subtree <- tree[c(1,3)]
  subranges <- ranges[c(1,3)]
  subpartition <- partition[c(1,3)]
  
  checkIdentical(as(subtree,"IRanges"), subranges)
  checkIdentical(subtree@partition, subpartition)
}

test_IntervalForest_length <- function() {
  ranges <- IRanges(c(1, 4, 9), c(5, 7, 10))
  tree <- IntervalForest(ranges, factor(c("a","a","b")))
  checkIdentical(length(tree), length(ranges))
}

test_IntervalForest_shift <- function() {
  ranges <- IRanges(c(1, 4, 9), c(5, 7, 10))
  partition <- Rle(factor(c("a","a","b")))
  tree <- IntervalForest(ranges, partition)
  tree <- shift(tree, 10)
  checkIdentical(as(tree, "IRanges"), shift(ranges, 10))
  checkIdentical(tree@partition, partition)
}