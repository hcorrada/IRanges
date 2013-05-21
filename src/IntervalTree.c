#include "IRanges.h"

#include "common.h"
#include "rbTree.h"
#include "localmem.h"

#include "ucsc_handlers.h"

#include <S.h> /* for Salloc() */

#define MAX(a, b)  (((a) > (b)) ? (a) : (b))

#define FIND_ALL    1
#define FIND_ANY    2
#define FIND_ARBITRARY 3

typedef struct _IntegerInterval {
  int start;
  int end;
} IntegerInterval;

typedef struct _IntegerIntervalNode {
  IntegerInterval interval;
  unsigned int index; /* identifies the interval */
  int maxEnd; /* maximum end value in children */
  unsigned int indexPosition;
} IntegerIntervalNode;

typedef struct _IntegerIntervalForest {
  struct rbTree **trees;
  int npartitions;
  int n;
  struct lm *lm;
} IntervalForest;


struct rbTree *_IntegerIntervalForest_getTree(IntervalForest *forest, int i) {
  return forest->trees[i];  
}

static int compare_interval(void *a, void *b) {
  IntegerInterval *ra = (IntegerInterval *)a;
  IntegerInterval *rb = (IntegerInterval *)b;
  if (rb->start < ra->start)
    return 1;
  /*
  if (ra->start < rb->start)
    return -1;
  if (rb->start == ra->start && ra->end != rb->end)
    return -1;
  return 0;
  */
  return -1;
}

static void _IntegerIntervalTree_free(SEXP r_tree) {
  struct rbTree *tree = (struct rbTree *)R_ExternalPtrAddr(r_tree);
  if (tree) {
    pushRHandlers();
    rbTreeFree(&tree);
    popRHandlers();
    R_ClearExternalPtr(r_tree);
  }
}


static void _IntegerIntervalForest_free(SEXP r_forest) {
  int npartitions;
  IntervalForest *forest = (IntervalForest *)R_ExternalPtrAddr(r_forest);
  struct rbTree *tree;
  int i;
  
  if (forest) {
    npartitions = forest->npartitions;
    pushRHandlers();
    for (i = 0; i<npartitions; i++) {
      tree = _IntegerIntervalForest_getTree(forest, i);
      freez(&tree);
    }
    lmCleanup(&forest->lm);
    popRHandlers();
    R_ClearExternalPtr(r_forest);
  }
}

static void _IntegerIntervalTree_add(struct rbTree *tree, int start, int end,
                       unsigned int index) {
  IntegerIntervalNode tmpInterval = { { start, end }, index, 0, index };
  IntegerIntervalNode *interval =
    lmCloneMem(tree->lm, &tmpInterval, sizeof(tmpInterval));
  rbTreeAdd(tree, interval);
}

static void _IntegerIndexedIntervalTree_add(struct rbTree *tree, int start, int end, unsigned int index, unsigned int indexPosition) {
  IntegerIntervalNode tmpInterval = { { start, end }, index, 0, indexPosition };
  IntegerIntervalNode *interval = lmCloneMem(tree->lm, &tmpInterval, sizeof(tmpInterval));
  rbTreeAdd(tree, interval);
}

/* need non-recursive implementation */
static void _IntegerIntervalNode_calc_max_end(struct rbTreeNode *node) {
  int maxEnd;
  maxEnd = ((IntegerInterval *)node->item)->end;
  if (node->left)
    _IntegerIntervalNode_calc_max_end(node->left);
  if (node->right)
    _IntegerIntervalNode_calc_max_end(node->right);
  if (node->left && node->right) {
    int childMax = MAX(((IntegerIntervalNode *)node->left->item)->maxEnd,
                       ((IntegerIntervalNode *)node->right->item)->maxEnd);
    maxEnd = MAX(childMax, maxEnd);
  } else if (node->left)
    maxEnd = MAX(((IntegerIntervalNode *)node->left->item)->maxEnd, maxEnd);
  else if (node->right)
    maxEnd = MAX(((IntegerIntervalNode *)node->right->item)->maxEnd, maxEnd);
  ((IntegerIntervalNode *)node->item)->maxEnd = maxEnd;
}

static void _IntegerIntervalTree_calc_max_end(struct rbTree *tree) {
  if (tree->root)
    _IntegerIntervalNode_calc_max_end(tree->root);
}

struct rbTree *_IntegerIntervalTree_new(void) {
  return rbTreeNew(compare_interval);
}

SEXP IntegerIntervalTree_new(SEXP r_ranges) {
  struct rbTree *tree = _IntegerIntervalTree_new();
  cachedIRanges cached_r_ranges = _cache_IRanges(r_ranges);
  int nranges = _get_cachedIRanges_length(&cached_r_ranges);
  int i, start, end;
  SEXP r_tree;

  pushRHandlers();
  for (i = 0; i < nranges; i++) {
    start = _get_cachedIRanges_elt_start(&cached_r_ranges, i);
    end = _get_cachedIRanges_elt_end(&cached_r_ranges, i);
    if (end >= start) /* only add non-empty ranges */
      _IntegerIntervalTree_add(tree, start, end, i+1);
  }
  popRHandlers();
  tree->n = nranges; /* kind of a hack - includes empty ranges */

  _IntegerIntervalTree_calc_max_end(tree);
  r_tree = R_MakeExternalPtr(tree, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(r_tree, _IntegerIntervalTree_free);

  return r_tree;
}

IntervalForest *_IntegerIntervalForest_new(int npartitions) {
  int i;
  struct lm *lm = lmInit(0);
  
  IntervalForest *forest;
  AllocVar(forest);
  
  forest->npartitions = npartitions;
  forest->lm = lm;
  forest->n = 0;
  
  lmAllocArray(forest->lm, forest->trees, npartitions * sizeof(struct rbTree *));
  for (i = 0; i<npartitions; i++) {
    struct rbTreeNode **stack = lmAlloc(forest->lm, 128 * sizeof(stack[0]));
    forest->trees[i] = rbTreeNewDetailed(compare_interval, forest->lm, stack);
  }
  return forest;
}

SEXP IntegerIntervalForest_new(SEXP r_ranges, SEXP r_partition, SEXP r_npartitions) {
  SEXP r_forest;
  
  cachedIRanges cached_r_ranges = _cache_IRanges(r_ranges);
  int nranges = _get_cachedIRanges_length(&cached_r_ranges);
  int i, start, end, cur_partition;
  int npartitions = *INTEGER(r_npartitions);
  int *partition = INTEGER(r_partition);
  
  IntervalForest *forest = _IntegerIntervalForest_new(npartitions);
  
  int tree_sizes[npartitions];
  for (i = 0; i < npartitions; i++) {
    tree_sizes[i] = 0;
  }
  struct rbTree *tree;

  pushRHandlers();
  for (i = 0; i < nranges; i++) {
    start = _get_cachedIRanges_elt_start(&cached_r_ranges, i);
    end = _get_cachedIRanges_elt_end(&cached_r_ranges, i);

    
    cur_partition = partition[i]-1;
    tree = _IntegerIntervalForest_getTree(forest, cur_partition);
    
    tree_sizes[cur_partition] += 1;    
    
    if (end >= start)
      _IntegerIntervalTree_add(tree, start, end, i+1);
  }
  popRHandlers();
  
  for (i = 0; i < npartitions; i++) {
    tree = _IntegerIntervalForest_getTree(forest, i);
    tree->n = tree_sizes[i];
    forest->n += tree_sizes[i];
    _IntegerIntervalTree_calc_max_end(tree);
  }
  
  r_forest = R_MakeExternalPtr(forest, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(r_forest, _IntegerIntervalForest_free);
  
  return r_forest;
}

/* If result_ints is NULL, returns a logical vector of with
 * length(r_ranges) elements whose values indicate whether or at
 * least one overlap was found.
 * If result_ints != NULL, adds every hit to result_ints and returns
 * an integer vector of with length(r_ranges) + 1 elements, where
 * the first length(r_ranges) contain a 0-based "partioning by end"
 * of overlap hit indices and the last element is the total number
 * of overlaps.
 * Running time: less than O(mk + mlogn) */
SEXP _IntegerIntervalTree_overlap(struct rbTree *tree, SEXP r_ranges,
                                  int find_type,
                                  struct slRef **result_ints)
{
  struct rbTreeNode *p = tree->root;
  struct slRef *active = NULL, *active_head = NULL;
  int /* stack height */ height = 0, /* query counter */ m;
  int /* result element */ *r_elt;
  Rboolean hit = FALSE; /* already hit current node? */
  SEXP result_inds;

  cachedIRanges cached_r_ranges = _cache_IRanges(r_ranges);
  int nranges = _get_cachedIRanges_length(&cached_r_ranges);
  if (find_type == FIND_ALL) {
    PROTECT(result_inds = allocVector(INTSXP, nranges + 1));
  } else if (find_type == FIND_ANY) {
    PROTECT(result_inds = allocVector(LGLSXP, nranges));
  } else if (find_type == FIND_ARBITRARY) {
    PROTECT(result_inds = allocVector(INTSXP, nranges));
  }
  memset(INTEGER(result_inds), 0, LENGTH(result_inds) * sizeof(int));
  
  if (!p) { /* tree is empty */
    UNPROTECT(1);
    return result_inds;
  }

  for (m = 0, r_elt = INTEGER(result_inds); m < nranges; m++, r_elt++) {
    int start = _get_cachedIRanges_elt_start(&cached_r_ranges, m);
    int end = _get_cachedIRanges_elt_end(&cached_r_ranges, m);
    if (end < start) { /* empty query range */
      if (find_type == FIND_ALL)
        *(r_elt+1) = *r_elt;
      continue;
    }
    int count = 0;
    /* add hits from previous query, if still valid */
    /* this trick lets us avoid restarting the search for every query */
    if (find_type == FIND_ALL) {
      struct slRef *prev = NULL;
      for (struct slRef *a = active_head; a != NULL;) {
        IntegerInterval *interval = a->val;
        if (interval->end < start) { /* never see this again */
          /* Rprintf("LINE %d -- goodbye: %d\n", __LINE__, ((IntegerIntervalNode*)interval)->index); */
          struct slRef *next = a->next;
          if (prev)
            prev->next = next;
          else
            active_head = next;
          freeMem(a);
          a = next;
        } else {
          if (interval->start > end) /* no more hits here */
            break;
          struct slRef *resultNode = slRefNew(interval);
          /* Rprintf("LINE %d -- p hit: %d\n", __LINE__, ((IntegerIntervalNode*)interval)->index); */
          slAddHead(result_ints, resultNode); /* owns Node */
          count++;
          prev = a;
          a = a->next;
        }
      }
      active = prev; /* active is the tail of the list (when it matters) */
    }
    while(1) {
      IntegerInterval *interval = (IntegerInterval *)p->item;
      /* is node on top of stack? */
      Rboolean visited = height && p == tree->stack[height-1];
      /* have to retry nodes on stack after query switch */
      /* Rprintf("LINE %d -- subject: %d,%d,%d / query: %d,%d, stack: %d\n", __LINE__,
              interval->start, interval->end,
              ((IntegerIntervalNode *)interval)->maxEnd, start, end, height); */

      /* in-order traversal of tree */

      /* go left if node not yet visited and max end satisfied */
      if(p->left && !visited && 
         ((IntegerIntervalNode *)p->left->item)->maxEnd >= start) {
        tree->stack[height++] = p;
        p = p->left;
        /*Rprintf("left\n");*/
      } else {
        /* consider current node if not already checked */
        if (interval->start <= end && interval->end >= start) {
          /* Rprintf("hit: %d\n", ((IntegerIntervalNode *)interval)->index); */
          if (find_type == FIND_ALL) {
            if (!hit) {
              struct slRef *resultNode = slRefNew(interval);
              slAddHead(result_ints, resultNode);
              resultNode = slRefNew(interval);
              if (active == NULL)
                active_head = resultNode;
              else
                active->next = resultNode;
              active = resultNode;
              count++;
            }
          } else if (find_type == FIND_ANY) {
            *r_elt = 1;
            break;
          } else if (find_type == FIND_ARBITRARY) {
            *r_elt = ((IntegerIntervalNode *)interval)->index;
            break;
          }
          hit = TRUE;
        }
        if (visited) { /* pop already visited node */ 
          height--;
          /*if (dirty) *//* retried this one */
          /*  dirty_level--;*/
        }

        /* go right if sensible */
        if (p->right && interval->start <= end &&
            ((IntegerIntervalNode *)p->right->item)->maxEnd >= start) {
          /* Rprintf("right\n"); */
          p = p->right;
          hit = FALSE;
        }
        else if (interval->start > end || !height) { /* no more hits */
          if (visited)
            height++; /* repush -- don't go left again */
          if (find_type == FIND_ALL)
            *(r_elt+1) = *r_elt + count;
          break;
        } else {
          p = tree->stack[height-1]; /* return to ancestor */
          hit = FALSE;
          /* Rprintf("up\n"); */
        }
      }
    }
  }
  UNPROTECT(1);
  return result_inds;
}

SEXP _IntegerIntervalForest_overlap(IntervalForest *forest, SEXP r_ranges, SEXP r_partition, int find_type, struct slRef **result_ints)
{
  int i, cur_partition = -1;
  struct rbTree *tree;
  struct rbTreeNode *p;
  struct slRef *active = NULL;
  struct slRef *active_head = NULL;
  /* stack height */
  int height = 0;
  
  /* already hit current node? */
  Rboolean hit = FALSE; 
  
  int /* query counter */ m;
  int /* result element */ *r_elt, /* partition element */ *p_elt;
  
  
  SEXP result_inds;
  
  cachedIRanges cached_r_ranges = _cache_IRanges(r_ranges);
  int nranges = _get_cachedIRanges_length(&cached_r_ranges);
  if (find_type == FIND_ALL) {
    PROTECT(result_inds = allocVector(INTSXP, nranges + 1));
  } else if (find_type == FIND_ANY) {
    PROTECT(result_inds = allocVector(LGLSXP, nranges));
  } else if (find_type == FIND_ARBITRARY) {
    PROTECT(result_inds = allocVector(INTSXP, nranges));
  }
  memset(INTEGER(result_inds), 0, LENGTH(result_inds) * sizeof(int));
  
  if (!forest->n) {
    UNPROTECT(1);
    return result_inds;
  }
  
  for (m = 0, r_elt = INTEGER(result_inds), p_elt = INTEGER(r_partition); m < nranges; m++, r_elt++, p_elt++) {
    int start = _get_cachedIRanges_elt_start(&cached_r_ranges, m);
    int end = _get_cachedIRanges_elt_end(&cached_r_ranges, m);
    if (end < start) { /* empty query range */
      if (find_type == FIND_ALL)
        *(r_elt+1) = *r_elt;
      continue;
    }
    
    if (*p_elt == NA_INTEGER) {
      /* query partition not in subject */
      if (find_type == FIND_ALL)
        *(r_elt+1) = *r_elt;
      continue;
    }
    if (cur_partition != *p_elt - 1) {
      /* new partition, reset everything */
      cur_partition = *p_elt - 1;
      tree = forest->trees[cur_partition];
      p = tree->root;
      active = NULL;
      active_head = NULL;
      height = 0;
      hit = FALSE;
    }
    
    /* Rprintf("LINE: %d -- cur_partition: %d, height: %d, hit: %d\n", __LINE__, cur_partition, height, hit); */
    if (!p) {
      /* this tree is empty */
      if (find_type == FIND_ALL)
        *(r_elt+1) = *r_elt;
      continue;
    }
    
    int count = 0;
    /* add hits from previous query, if still valid */
    /* this trick lets us avoid restarting the search for every query */
    if (find_type == FIND_ALL) {
      struct slRef *prev = NULL;
      for (struct slRef *a = active_head; a != NULL;) {
        IntegerInterval *interval = a->val;
        if (interval->end < start) { /* never see this again */
          /* Rprintf("LINE: %d -- goodbye: %d\n", __LINE__, ((IntegerIntervalNode*)interval)->index); */
          struct slRef *next = a->next;
          if (prev)
            prev->next = next;
          else
            active_head = next;
          freeMem(a);
          a = next;
        } else {
          if (interval->start > end) /* no more hits here */
            break;
          struct slRef *resultNode = slRefNew(interval);
          /* Rprintf("LINE: %d -- p hit: %d\n", __LINE__, ((IntegerIntervalNode*)interval)->index); */
          slAddHead(result_ints, resultNode); /* owns Node */
          count++;
          prev = a;
          a = a->next;
        }
      }
      active = prev; /* active is the tail of the list (when it matters) */
    }
    while(1) {
      IntegerInterval *interval = (IntegerInterval *)p->item;
      /* is node on top of stack? */
      Rboolean visited = height && p == tree->stack[height-1];
      /* have to retry nodes on stack after query switch */
      /* Rprintf("LINE %d -- subject: %d,%d,%d / query: %d,%d, stack: %d\n", __LINE__,
              interval->start, interval->end,
              ((IntegerIntervalNode *)interval)->maxEnd, start, end, height); */

      /* in-order traversal of tree */

      /* go left if node not yet visited and max end satisfied */
      if(p->left && !visited && 
         ((IntegerIntervalNode *)p->left->item)->maxEnd >= start) {
        tree->stack[height++] = p;
        p = p->left;
        /* Rprintf("LINE %d -- left\n", __LINE__); */
      } else {
         /*consider current node if not already checked */
        if (interval->start <= end && interval->end >= start) {
           /* Rprintf("LINE %d -- hit: %d\n", __LINE__, ((IntegerIntervalNode *)interval)->index); */
          if (find_type == FIND_ALL) {
            if (!hit) {
              struct slRef *resultNode = slRefNew(interval);
              slAddHead(result_ints, resultNode);
              resultNode = slRefNew(interval);
              if (active == NULL)
                active_head = resultNode;
              else
                active->next = resultNode;
              active = resultNode;
              count++;
            }
          } else if (find_type == FIND_ANY) {
            *r_elt = 1;
            break;
          } else if (find_type == FIND_ARBITRARY) {
            *r_elt = ((IntegerIntervalNode *)interval)->index;
            break;
          }
          hit = TRUE;
        }
        if (visited) { /* pop already visited node */ 
          height--;
          /*if (dirty) *//* retried this one */
          /*  dirty_level--;*/
        }

        /* go right if sensible */
        if (p->right && interval->start <= end &&
            ((IntegerIntervalNode *)p->right->item)->maxEnd >= start) {
           /* Rprintf("LINE %d -- right\n", __LINE__); */
          p = p->right;
          hit = FALSE;
        }
        else if (interval->start > end || !height) { /* no more hits */
          if (visited)
            height++; /* repush -- don't go left again */
          if (find_type == FIND_ALL)
            *(r_elt+1) = *r_elt + count;
          break;
        } else {
          p = tree->stack[height-1]; /* return to ancestor */
          hit = FALSE;
           /*Rprintf("LINE %d -- up\n", __LINE__); */
        }
      }
    }
    /*Rprintf("LINE %d -- height: %d, hit: %d\n", __LINE__, height, hit);*/
  }
  
  UNPROTECT(1);
  return result_inds;
}

SEXP IntegerIntervalTree_overlap_any(SEXP r_tree, SEXP r_ranges,
                                     SEXP r_order) {
  int i, *left, *right, *o_elt, nranges = _get_IRanges_length(r_ranges);
  struct rbTree *tree = R_ExternalPtrAddr(r_tree);
  pushRHandlers();
  SEXP r_unordered =
    _IntegerIntervalTree_overlap(tree, r_ranges, FIND_ANY, NULL);
  popRHandlers();
  PROTECT(r_unordered);
  SEXP r_ordered = allocVector(LGLSXP, nranges);
  left = INTEGER(r_ordered);
  for (i = 0, right = INTEGER(r_unordered), o_elt = INTEGER(r_order);
       i < nranges; i++, right++, o_elt++) {
    left[*o_elt - 1] = *right; 
  }
  UNPROTECT(1);
  return r_ordered;
}

SEXP IntegerIntervalForest_overlap_any(SEXP r_forest, SEXP r_ranges, SEXP r_partition, SEXP r_order)
{
  int i, *left, *right, *o_elt, nranges = _get_IRanges_length(r_ranges);
  IntervalForest *forest = R_ExternalPtrAddr(r_forest);
  pushRHandlers();
  SEXP r_unordered =
    _IntegerIntervalForest_overlap(forest, r_ranges, r_partition, FIND_ANY, NULL);
  popRHandlers();
  PROTECT(r_unordered);
  SEXP r_ordered = allocVector(LGLSXP, nranges);
  left = INTEGER(r_ordered);
  for (i = 0, right = INTEGER(r_unordered), o_elt = INTEGER(r_order);
       i < nranges; i++, right++, o_elt++) {
    left[*o_elt - 1] = *right; 
  }
  UNPROTECT(1);
  return r_ordered;
}

SEXP IntegerIntervalTree_overlap_arbitrary(SEXP r_tree, SEXP r_ranges,
                                           SEXP r_order) {
  int i, *left, *right, *o_elt, nranges = _get_IRanges_length(r_ranges);
  struct rbTree *tree = R_ExternalPtrAddr(r_tree);
  pushRHandlers();
  SEXP r_unordered =
    _IntegerIntervalTree_overlap(tree, r_ranges, FIND_ARBITRARY, NULL);
  popRHandlers();
  PROTECT(r_unordered);
  SEXP r_ordered = allocVector(INTSXP, nranges);
  left = INTEGER(r_ordered);
  for (i = 0, right = INTEGER(r_unordered), o_elt = INTEGER(r_order);
       i < nranges; i++, right++, o_elt++) {
    left[*o_elt - 1] = *right > 0 ? *right : NA_INTEGER;
  }
  UNPROTECT(1);
  return r_ordered;
}

SEXP IntegerIntervalForest_overlap_arbitrary(SEXP r_forest, SEXP r_ranges, SEXP r_partition, SEXP r_order)
{
  int i, *left, *right, *o_elt, nranges = _get_IRanges_length(r_ranges);
  IntervalForest *forest = R_ExternalPtrAddr(r_forest);
  pushRHandlers();
  SEXP r_unordered =
    _IntegerIntervalForest_overlap(forest, r_ranges, r_partition, FIND_ARBITRARY, NULL);
  popRHandlers();
  PROTECT(r_unordered);
  SEXP r_ordered = allocVector(INTSXP, nranges);
  left = INTEGER(r_ordered);
  for (i = 0, right = INTEGER(r_unordered), o_elt = INTEGER(r_order);
       i < nranges; i++, right++, o_elt++) {
    left[*o_elt - 1] = *right > 0 ? *right : NA_INTEGER;
  }
  UNPROTECT(1);
  return r_ordered;
}
SEXP IntegerIntervalTree_overlap_first(SEXP r_tree, SEXP r_ranges,
                                       SEXP r_order)
{
  struct rbTree *tree = R_ExternalPtrAddr(r_tree);
  struct slRef *results = NULL, *result;
  SEXP r_query_start, r_results;
  int i, j, index, nhits, nranges = _get_IRanges_length(r_ranges);
  int *left, *right, *r_vector, *r_elt, *o_elt;

  pushRHandlers();
  r_query_start =
    _IntegerIntervalTree_overlap(tree, r_ranges, FIND_ALL, &results);
  PROTECT(r_query_start);
  nhits = INTEGER(r_query_start)[nranges];
  slReverse(&results);

  PROTECT(r_results = allocVector(INTSXP, nranges));
  for (i = 0, r_elt = INTEGER(r_results); i < nranges; i++, r_elt++)
    *r_elt = NA_INTEGER;
  result = results;
  r_vector = INTEGER(r_results);
  for (i = 0, o_elt = INTEGER(r_order),
       left = INTEGER(r_query_start), right = INTEGER(r_query_start) + 1;
       i < nranges; i++, o_elt++, left++, right++) {
    r_elt = r_vector + (*o_elt - 1);
    for (j = *left; j < *right; j++) {
      index = ((IntegerIntervalNode *)result->val)->indexPosition;
      if (*r_elt == NA_INTEGER || (*r_elt > index))
        *r_elt = index;
      result = result->next;
    }
  }

  slFreeList(&results);
  popRHandlers();
  
  UNPROTECT(2);

  return r_results;
}

SEXP IntegerIntervalForest_overlap_first(SEXP r_forest, SEXP r_ranges, SEXP r_partition, SEXP r_order)
{
  IntervalForest *forest = R_ExternalPtrAddr(r_forest);
  struct slRef *results = NULL, *result;
  SEXP r_query_start, r_results;
  int i, j, index, nhits, nranges = _get_IRanges_length(r_ranges);
  int *left, *right, *r_vector, *r_elt, *o_elt;
  
  pushRHandlers();
  r_query_start = _IntegerIntervalForest_overlap(forest, r_ranges, r_partition, FIND_ALL, &results);
  PROTECT(r_query_start);
  nhits = INTEGER(r_query_start)[nranges];
  slReverse(&results);
  
  PROTECT(r_results = allocVector(INTSXP, nranges));
  for (i = 0, r_elt = INTEGER(r_results); i < nranges; i++, r_elt++)
    *r_elt = NA_INTEGER;
  result = results;
  r_vector = INTEGER(r_results);
  for (i = 0, o_elt = INTEGER(r_order), left = INTEGER(r_query_start), right=INTEGER(r_query_start) + 1;
       i < nranges; i++, o_elt++, left++, right++) {
         r_elt = r_vector + (*o_elt - 1);
         
         for (j = *left; j < *right; j++) {
           index = ((IntegerIntervalNode *)result->val)->index;
           if (*r_elt == NA_INTEGER || (*r_elt > index)) {
              *r_elt = index;
           }
          result = result->next;
        }
  }
  
  slFreeList(&results);
  popRHandlers();
  UNPROTECT(2);
  
  return r_results;
}

SEXP IntegerIntervalTree_overlap_last(SEXP r_tree, SEXP r_ranges,
                                      SEXP r_order)
{
  struct rbTree *tree = R_ExternalPtrAddr(r_tree);
  struct slRef *results = NULL, *result;
  SEXP r_query_start, r_results;
  int i, j, index, nhits, nranges = _get_IRanges_length(r_ranges);
  int *left, *right, *r_vector, *r_elt, *o_elt;

  pushRHandlers();
  r_query_start =
    _IntegerIntervalTree_overlap(tree, r_ranges, FIND_ALL, &results);
  PROTECT(r_query_start);
  nhits = INTEGER(r_query_start)[nranges];
  slReverse(&results);

  PROTECT(r_results = allocVector(INTSXP, nranges));
  for (i = 0, r_elt = INTEGER(r_results); i < nranges; i++, r_elt++)
    *r_elt = NA_INTEGER;
  result = results;
  r_vector = INTEGER(r_results);
  for (i = 0, o_elt = INTEGER(r_order),
       left = INTEGER(r_query_start), right = INTEGER(r_query_start) + 1;
       i < nranges; i++, o_elt++, left++, right++) {
    r_elt = r_vector + (*o_elt - 1);
    for (j = *left; j < *right; j++) {
      index = ((IntegerIntervalNode *)result->val)->indexPosition;
      if (*r_elt == NA_INTEGER || (*r_elt < index))
        *r_elt = index;
      result = result->next;
    }
  }

  slFreeList(&results);
  popRHandlers();
  
  UNPROTECT(2);

  return r_results;
}

SEXP IntegerIntervalForest_overlap_last(SEXP r_forest, SEXP r_ranges, SEXP r_partition, SEXP r_order)
{
  IntervalForest *forest = R_ExternalPtrAddr(r_forest);
  struct slRef *results = NULL, *result;
  SEXP r_query_start, r_results;
  int i, j, index, nhits, nranges = _get_IRanges_length(r_ranges);
  int *left, *right, *r_vector, *r_elt, *o_elt;

  pushRHandlers();
  r_query_start =
    _IntegerIntervalForest_overlap(forest, r_ranges, r_partition, FIND_ALL, &results);
  PROTECT(r_query_start);
  nhits = INTEGER(r_query_start)[nranges];
  slReverse(&results);

  PROTECT(r_results = allocVector(INTSXP, nranges));
  for (i = 0, r_elt = INTEGER(r_results); i < nranges; i++, r_elt++)
    *r_elt = NA_INTEGER;
  result = results;
  r_vector = INTEGER(r_results);
  for (i = 0, o_elt = INTEGER(r_order),
       left = INTEGER(r_query_start), right = INTEGER(r_query_start) + 1;
       i < nranges; i++, o_elt++, left++, right++) {
    r_elt = r_vector + (*o_elt - 1);
    for (j = *left; j < *right; j++) {
      index = ((IntegerIntervalNode *)result->val)->index;
      if (*r_elt == NA_INTEGER || (*r_elt < index))
        *r_elt = index;
      result = result->next;
    }
  }

  slFreeList(&results);
  popRHandlers();
  
  UNPROTECT(2);

  return r_results;
}

SEXP IntegerIntervalTree_overlap_all(SEXP r_tree, SEXP r_ranges, SEXP r_order)
{
  struct rbTree *tree = R_ExternalPtrAddr(r_tree);
  struct slRef *results = NULL, *result;
  SEXP r_query_start, r_results, r_query_hits, r_subject_hits;
  int i, j, nhits, nranges = _get_IRanges_length(r_ranges);
  int *left, *right, *r_elt, *o_elt;

  pushRHandlers();
  r_query_start =
    _IntegerIntervalTree_overlap(tree, r_ranges, FIND_ALL, &results);
  PROTECT(r_query_start);
  nhits = INTEGER(r_query_start)[nranges];
  slReverse(&results);

  int *r_query_col = (int *) R_alloc((long) nhits, sizeof(int));
  r_elt = r_query_col;
  for (i = 1, o_elt = INTEGER(r_order),
       left = INTEGER(r_query_start), right = INTEGER(r_query_start) + 1;
       i < LENGTH(r_query_start); i++, o_elt++, left++, right++) {
    for (j = *left; j < *right; j++) {
      *r_elt = *o_elt;
      r_elt++;
    }
  }

  int *r_subject_col = (int *) R_alloc((long) nhits, sizeof(int));
  for (result = results, r_elt = r_subject_col; result != NULL;
       result = result->next, r_elt++)
    *r_elt = ((IntegerIntervalNode *)result->val)->indexPosition;

  int *row = (int *) R_alloc((long) nhits, sizeof(int));
  _get_order_of_int_pairs(r_query_col, r_subject_col, nhits, 0, row, 0);

  PROTECT(r_results = NEW_OBJECT(MAKE_CLASS("Hits")));

  r_query_hits = NEW_INTEGER(nhits);
  SET_SLOT(r_results, install("queryHits"), r_query_hits);
  r_subject_hits = NEW_INTEGER(nhits);
  SET_SLOT(r_results, install("subjectHits"), r_subject_hits);

  for (i = 0, left = INTEGER(r_query_hits), right = INTEGER(r_subject_hits),
       o_elt = row; i < nhits; i++, left++, right++, o_elt++) {
    *left = r_query_col[*o_elt];
    *right = r_subject_col[*o_elt];
  }

  SET_SLOT(r_results, install("queryLength"), ScalarInteger(nranges));
  SET_SLOT(r_results, install("subjectLength"), ScalarInteger(tree->n));

  slFreeList(&results);
  popRHandlers();

  UNPROTECT(2);

  return r_results;
}

SEXP IntegerIntervalForest_overlap_all(SEXP r_forest, SEXP r_ranges, SEXP r_partition, SEXP r_order)
{
  IntervalForest *forest = R_ExternalPtrAddr(r_forest);
  struct slRef *results = NULL, *result;
  SEXP r_query_start, r_results, r_query_hits, r_subject_hits;
  int i, j, nhits, nranges = _get_IRanges_length(r_ranges);
  int *left, *right, *r_elt, *o_elt;

  pushRHandlers();
  r_query_start =
    _IntegerIntervalForest_overlap(forest, r_ranges, r_partition, FIND_ALL, &results);
  PROTECT(r_query_start);
  nhits = INTEGER(r_query_start)[nranges];
  slReverse(&results);

  int *r_query_col = (int *) R_alloc((long) nhits, sizeof(int));
  r_elt = r_query_col;
  for (i = 1, o_elt = INTEGER(r_order),
       left = INTEGER(r_query_start), right = INTEGER(r_query_start) + 1;
       i < LENGTH(r_query_start); i++, o_elt++, left++, right++) {
    for (j = *left; j < *right; j++) {
      *r_elt = *o_elt;
      r_elt++;
    }
  }

  int *r_subject_col = (int *) R_alloc((long) nhits, sizeof(int));
  for (result = results, r_elt = r_subject_col; result != NULL;
       result = result->next, r_elt++)
    *r_elt = ((IntegerIntervalNode *)result->val)->index;

  int *row = (int *) R_alloc((long) nhits, sizeof(int));
  _get_order_of_int_pairs(r_query_col, r_subject_col, nhits, 0, row, 0);

  PROTECT(r_results = NEW_OBJECT(MAKE_CLASS("Hits")));

  r_query_hits = NEW_INTEGER(nhits);
  SET_SLOT(r_results, install("queryHits"), r_query_hits);
  r_subject_hits = NEW_INTEGER(nhits);
  SET_SLOT(r_results, install("subjectHits"), r_subject_hits);

  for (i = 0, left = INTEGER(r_query_hits), right = INTEGER(r_subject_hits),
       o_elt = row; i < nhits; i++, left++, right++, o_elt++) {
    *left = r_query_col[*o_elt];
    *right = r_subject_col[*o_elt];
  }

  SET_SLOT(r_results, install("queryLength"), ScalarInteger(nranges));
  SET_SLOT(r_results, install("subjectLength"), ScalarInteger(forest->n));

  slFreeList(&results);
  popRHandlers();

  UNPROTECT(2);

  return r_results;
}

// finds the node to the right of each element of sorted query
/*
struct rbTreeNode**
_IntegerIntervalTree_nearest(struct rbTree *tree, int *query, int len,
                             Rboolean self)
{
  struct rbTreeNode *p, *nextP;
  int i = 0;
  rbTreeNode **result = R_alloc(len, sizeof(struct rbTreeNode *));
  for (i = 0, p = tree->root; p != NULL && i < len; p = nextP) {
    IntegerInterval *subject = (IntegerInterval *)p->item;
    Rboolean visited = height && p == tree->stack[height-1];  
    int sep = subject->start - query[i];
    nextP = NULL;
    if(sep > 0 && !visited) { // can go left
      tree->stack[height++] = p;
      nextP = p->left;
    } else {
      if (visited)
        height--; // pop handled node if on stack
      if(sep <= 0)
        nextP = p->right; 
    }
    if (!nextP) {
      nextP = height ? tree->stack[height-1] : p;
      IntegerInterval *other = (IntegerInterval *)nextP->item;
      result[i++] = other->start;
    }
  }
  
  return result_ind;
}
*/

/*
SEXP IntegerIntervalTree_nearest(SEXP r_tree, SEXP r_ranges, SEXP r_self) {
  struct rbTree *tree = R_ExternalPtrAddr(r_tree);
  cachedIRanges cached_r_ranges = _cache_IRanges(r_ranges);
  int nranges = _get_cachedIRanges_length(&cached_r_ranges);


  int *end = R_alloc(nranges, sizeof(int));
    
  for (i = 0; i < nranges; i++)
    end[i] = _get_cachedIRanges_elt_end(&cached_r_ranges, i);
  
}
*/

/* Traverses the tree, pulling out the intervals, in order */
IntegerInterval **_IntegerIntervalTree_intervals(struct rbTree *tree) {
  struct rbTreeNode *p = tree->root;
  int height = 0;
  IntegerInterval **intervals = Salloc(tree->n, IntegerInterval *);

  if (tree->n && p)
    while(1) {
      /* is node on top of stack? */
      Rboolean visited = height && p == tree->stack[height-1];
      /* first, check for overlap */
      if (!visited && p->left) {
        /* push current node onto stack */
        tree->stack[height++] = p;
        /* go left */
        p = p->left;
      } else {       /* can't go left, handle this node */
        intervals[((IntegerIntervalNode *)p->item)->index - 1] =
          (IntegerInterval *)p->item;
        if (visited)
          height--; /* pop handled node if on stack */
        if (p->right) /* go right if possible */
          p = p->right;
        else if (height) /* more on stack */
          p = tree->stack[height-1];
        else
          break; /* nothing left on stack, we're finished */
      }
    }

  return intervals;
}

/* Traverses the forest, pulling out the intervals, in order */
/* this should be cleaned up as it duplicates a lot of code above */
IntegerInterval **_IntegerIntervalForest_intervals(IntervalForest *forest) {
  struct rbTree *tree;
  struct rbTreeNode *p;
  
  int height, cur_partition;
  IntegerInterval **intervals = Salloc(forest->n, IntegerInterval *);
  
  for (cur_partition = 0; cur_partition < forest->npartitions; cur_partition++) {
    tree = _IntegerIntervalForest_getTree(forest, cur_partition);
    p = tree->root;
    height = 0;
    
    if (tree->n && p) {
      while (1) {
        /* is node on top of stack? */
        Rboolean visited = height && p == tree->stack[height=1];
        /* first, check for overlap */
        if (!visited && p->left) {
          /* push current node onto stack */
          tree->stack[height++] = p;
          /* go left */
          p = p->left;
        } else { /* can't go left, handle this node */
          intervals[((IntegerIntervalNode *)p->item)->index - 1] = (IntegerInterval *)p->item;
          if (visited) {
            height--; /* pop handled node if on stack */
          }
          if (p->right) { /* for right if possible */
            p = p->right;
          }
          else if (height) { /* more on stack */
            p = tree->stack[height-1];
          } else {
            break; /* nothing left on stack, we're finished */
          }
        }
      }
    }
  }
  return intervals;
}


SEXP IntegerIntervalTree_asIRanges(SEXP r_tree) {
  SEXP r_start, r_width, r_ranges;
  struct rbTree *tree = R_ExternalPtrAddr(r_tree);
  pushRHandlers();
  IntegerInterval **intervals = _IntegerIntervalTree_intervals(tree);
  popRHandlers();
  int i, *s_elt, *w_elt;

  /* Rprintf("tree size %d\n", tree->n); */
  
  PROTECT(r_start = allocVector(INTSXP, tree->n));
  PROTECT(r_width = allocVector(INTSXP, tree->n));

  for(i = 0, s_elt = INTEGER(r_start), w_elt = INTEGER(r_width); i < tree->n;
      i++, s_elt++, w_elt++) {
    if (intervals[i]) {
      *s_elt = intervals[i]->start;
      *w_elt = intervals[i]->end - intervals[i]->start + 1;
    } else {
      *s_elt = 1;
      *w_elt = 0;
    }
  }

  r_ranges = _new_IRanges("IRanges", r_start, r_width, R_NilValue);

  UNPROTECT(2);
  return r_ranges;
}

SEXP IntegerIntervalTree_start(SEXP r_tree) {
  SEXP r_start;
  struct rbTree *tree = R_ExternalPtrAddr(r_tree);
  pushRHandlers();
  IntegerInterval **intervals = _IntegerIntervalTree_intervals(tree);
  popRHandlers();
  int i, *r_elt;

  r_start = allocVector(INTSXP, tree->n);

  for(i = 0, r_elt = INTEGER(r_start); i < tree->n; i++, r_elt++)
    *r_elt = intervals[i] ? intervals[i]->start : 1;

  return r_start;
}

SEXP IntegerIntervalForest_start(SEXP r_forest) {
  SEXP r_start;
  IntervalForest *forest = R_ExternalPtrAddr(r_forest);
  pushRHandlers();
  IntegerInterval **intervals = _IntegerIntervalForest_intervals(forest);
  popRHandlers();
  int i, *r_elt;
  
  r_start = allocVector(INTSXP, forest->n);  
    
  for (i = 0, r_elt = INTEGER(r_start); i < forest->n; i++, r_elt++) {
      *r_elt = intervals[i] ? intervals[i]->start : 1;
  }
  return r_start;
}

SEXP IntegerIntervalTree_end(SEXP r_tree) {
  SEXP r_end;
  struct rbTree *tree = R_ExternalPtrAddr(r_tree);
  pushRHandlers();
  IntegerInterval **intervals = _IntegerIntervalTree_intervals(tree);
  popRHandlers();
  int i, *r_elt;

  r_end = allocVector(INTSXP, tree->n);

  for(i = 0, r_elt = INTEGER(r_end); i < tree->n; i++, r_elt++)
    *r_elt = intervals[i] ? intervals[i]->end : 0;

  return r_end;
}

SEXP IntegerIntervalForest_end(SEXP r_forest) {
  SEXP r_end;
  IntervalForest *forest = R_ExternalPtrAddr(r_forest);
  pushRHandlers();
  IntegerInterval **intervals = _IntegerIntervalForest_intervals(forest);
  popRHandlers();
  int i, *r_elt;
  
  r_end = allocVector(INTSXP, forest->n);
  
  for (i = 0, r_elt = INTEGER(r_end); i < forest->n; i++, r_elt++) {
    *r_elt = intervals[i] ? intervals[i]->end : 0;
  }
  
  return r_end;
  return ScalarInteger(12345);  
}

SEXP IntegerIntervalTree_length(SEXP r_tree) {
  struct rbTree *tree = R_ExternalPtrAddr(r_tree);
  return ScalarInteger(tree->n);
}

SEXP IntegerIntervalForest_length(SEXP r_forest) {
  IntervalForest *forest = R_ExternalPtrAddr(r_forest);
  return ScalarInteger(forest->n);  
}

static void _IntegerIntervalTreeNode_dump(void *item, FILE *file) {
  IntegerInterval *node = (IntegerInterval *)item;
  fprintf(file, "index %d -- %d:%d / %d", ((IntegerIntervalNode *)node)->index, node->start, node->end,
          ((IntegerIntervalNode *)node)->maxEnd);
}

SEXP IntegerIntervalTree_dump(SEXP r_tree) {
  struct rbTree *tree = R_ExternalPtrAddr(r_tree);
  pushRHandlers();
  rbTreeDump(tree, stdout, _IntegerIntervalTreeNode_dump);
  popRHandlers();
  return R_NilValue;
}

SEXP IntegerIntervalForest_dump(SEXP r_forest) {
  IntervalForest *forest = R_ExternalPtrAddr(r_forest);
  for (int i = 0; i < forest->npartitions; i++) {
    struct rbTree *tree = forest->trees[i];
    fprintf(stdout, "Partition: %d length: %d\n", i, tree->n);
    pushRHandlers();
    rbTreeDump(tree, stdout, _IntegerIntervalTreeNode_dump);
    popRHandlers();
    fprintf(stdout, "\n");
  }
  return R_NilValue;
}

