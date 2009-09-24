/****************************************************************************
 *                   Low-level utilities for copying data                   *
 *                from a vector to a vector of the same type                *
 *                ------------------------------------------                *
 *                                                                          *
 * NOTE: Only "raw" (RAWSXP), "logical" (LGLSXP), "integer" (INTSXP),       *
 *       "double" (REALSXP) and "complex" (CPLXSXP) vectors are supported   *
 *       at the moment.                                                     *
 ****************************************************************************/
#include "IRanges.h"


/****************************************************************************
 * A. memcpy()-BASED COPY (NO RECYCLING)
 * =====================================
 */

void _vector_memcpy(SEXP out, int out_offset, SEXP in, int in_offset, int nelt)
{
	if (out_offset < 0 || out_offset + nelt > LENGTH(out)
	 || in_offset < 0 || in_offset + nelt > LENGTH(in))
		error("subscripts out of bounds");
	switch (TYPEOF(out)) {
	case RAWSXP:
		memcpy(RAW(out) + out_offset, RAW(in) + in_offset,
			nelt * sizeof(Rbyte));
		break;
	case LGLSXP:
	case INTSXP:
		memcpy(INTEGER(out) + out_offset, INTEGER(in) + in_offset,
			nelt * sizeof(int));
		break;
	case REALSXP:
		memcpy(REAL(out) + out_offset, REAL(in) + in_offset,
			nelt * sizeof(double));
		break;
	case CPLXSXP:
		memcpy(COMPLEX(out) + out_offset, COMPLEX(in) + in_offset,
			nelt * sizeof(Rcomplex));
		break;
	default:
		error("IRanges internal error in _vector_memcpy(): "
		      "%s type not supported", type2str(TYPEOF(out)));
	}
	return;
}


/****************************************************************************
 * B. COPY WITH RECYCLING
 * ======================
 *
 * The functions in section B. implement copy with recycling. The user can
 * choose between 2 interfaces for specifying elements in the 'in' or 'out'
 * vectors:
 *
 *   1. The "offset/nelt" interface: the elements to access are specified via
 * 2 integers: 'offset' (the 0-based position of the first element to access)
 * and 'nelt' (the number of elements to access, all immediately following
 * the first element to access).
 *
 *   2. The "subset" interface: the elements to access are specified by an
 * integer vector containing their 1-based positions in the 'in' or 'out'
 * vectors.
 *
 * The "subset" interface is intended to be used by the subsetting
 * operator [ defined at the R level for SharedVector objects.
 * Implementing this interface requires to pay some special attention to
 * the following important properties of the subsetting operator [ in R.
 * If x is a vector and i an integer vector of length n with the following
 * properties:
 *   a) i contains no NA values,
 *   b) i can be used to subset x without being "out of bounds" (i.e all
 *      values in i are >= 1 and <= length(x)),
 * then we have the following properties:
 *   1) READING from x: y <- x[i] produces a vector, of the same type than x,
 *      but of the same length than i (length(y) == n).
 *   2) READING from then WRITING to x: x[i] <- x[i] (short for y <- x[i];
 *      x[i] <- y) doesn't modify the values in x.
 *   3) WRITING to then READING from x: if z is a vector of length n and of
 *      the same type than x, then doing x[i] <- z; y <- x[i] guarantees that
 *      y is identical to z only when i contains no repeated value!
 *
 * Functions in this file that implement the "subset" interface adhere to the
 * above properties.
 */


/*
 * RECYCLING: Cyclic writing to 'out'.
 * INTERFACE: "offset/nelt".
 * In addition, "raw" vectors support fast on-the-fly translation via the
 * 'lkup' table.
 * Reverts the order of the copied elements if 'reverse' is != 0.
 */
void _vector_copy_from_offset(SEXP out, SEXP in, int in_offset, int nelt,
		SEXP lkup, int reverse)
{
	int i1, i2;
	void (*fun)(int, int, char *, size_t, const char *, size_t, size_t);

	i1 = in_offset;
	i2 = in_offset + nelt - 1;
	fun = reverse ? _IRanges_reverse_memcpy_from_i1i2
		      : _IRanges_memcpy_from_i1i2;
	switch (TYPEOF(out)) {
	case RAWSXP:
		if (lkup == R_NilValue) {
			fun(i1, i2,
				(char *) RAW(out), LENGTH(out),
				(char *) RAW(in), LENGTH(in), sizeof(Rbyte));
		} else {
			if (reverse)
				_IRanges_reverse_charcpy_from_i1i2_with_lkup(
					i1, i2,
					(char *) RAW(out), LENGTH(out),
					(char *) RAW(in), LENGTH(in),
					INTEGER(lkup), LENGTH(lkup));
			else
				_IRanges_charcpy_from_i1i2_with_lkup(i1, i2,
					(char *) RAW(out), LENGTH(out),
					(char *) RAW(in), LENGTH(in),
					INTEGER(lkup), LENGTH(lkup));
		}
		break;
	case LGLSXP:
	case INTSXP:
		fun(i1, i2,
			(char *) INTEGER(out), LENGTH(out),
			(char *) INTEGER(in), LENGTH(in), sizeof(int));
		break;
	case REALSXP:
		fun(i1, i2,
			(char *) REAL(out), LENGTH(out),
			(char *) REAL(in), LENGTH(in), sizeof(double));
		break;
	case CPLXSXP:
		fun(i1, i2,
			(char *) COMPLEX(out), LENGTH(out),
			(char *) COMPLEX(in), LENGTH(in), sizeof(Rcomplex));
		break;
	default:
		error("IRanges internal error in _vector_copy_from_offset(): "
		      "%s type not supported", type2str(TYPEOF(out)));
	}
	return;
}

/*
 * RECYCLING: Cyclic writing to 'out'.
 * INTERFACE: "subset".
 * In addition, "raw" vectors support fast on-the-fly translation via the
 * 'lkup' table.
 */
void _vector_copy_from_subset(SEXP out, SEXP in, SEXP subset, SEXP lkup)
{
	switch (TYPEOF(out)) {
	case RAWSXP:
		if (lkup == R_NilValue)
			_IRanges_memcpy_from_subset(
				INTEGER(subset), LENGTH(subset),
				(char *) RAW(out), LENGTH(out),
				(char *) RAW(in), LENGTH(in), sizeof(Rbyte));
		else
			_IRanges_charcpy_from_subset_with_lkup(
				INTEGER(subset), LENGTH(subset),
				(char *) RAW(out), LENGTH(out),
				(char *) RAW(in), LENGTH(in), 
				INTEGER(lkup), LENGTH(lkup));
		break;
	case LGLSXP:
	case INTSXP:
		_IRanges_memcpy_from_subset(INTEGER(subset), LENGTH(subset),
			(char *) INTEGER(out), LENGTH(out),
			(char *) INTEGER(in), LENGTH(in), sizeof(int));
		break;
	case REALSXP:
		_IRanges_memcpy_from_subset(INTEGER(subset), LENGTH(subset),
			(char *) REAL(out), LENGTH(out),
			(char *) REAL(in), LENGTH(in), sizeof(double));
		break;
	case CPLXSXP:
		_IRanges_memcpy_from_subset(INTEGER(subset), LENGTH(subset),
			(char *) COMPLEX(out), LENGTH(out),
			(char *) COMPLEX(in), LENGTH(in), sizeof(Rcomplex));
		break;
	default:
		error("IRanges internal error in _vector_copy_from_subset(): "
		      "%s type not supported", type2str(TYPEOF(out)));
	}
	return;
}

/*
 * RECYCLING: Cyclic reading from 'in'.
 * INTERFACE: "offset/nelt".
 *
 * In addition, "raw" vectors support fast on-the-fly translation via the
 * 'lkup' table.
 */
void _vector_copy_to_offset(SEXP out, SEXP in, int out_offset, int nelt,
		SEXP lkup)
{
	int i1, i2;

	i1 = out_offset;
	i2 = out_offset + nelt - 1;
	switch (TYPEOF(out)) {
	case RAWSXP:
		if (lkup == R_NilValue) {
			_IRanges_memcpy_to_i1i2(i1, i2,
				(char *) RAW(out), LENGTH(out),
				(char *) RAW(in), LENGTH(in), sizeof(Rbyte));
		} else {
			_IRanges_charcpy_to_i1i2_with_lkup(i1, i2,
				(char *) RAW(out), LENGTH(out),
				(char *) RAW(in), LENGTH(in),
				INTEGER(lkup), LENGTH(lkup));
		}
		break;
	case LGLSXP:
	case INTSXP:
		_IRanges_memcpy_to_i1i2(i1, i2,
			(char *) INTEGER(out), LENGTH(out),
			(char *) INTEGER(in), LENGTH(in), sizeof(int));
		break;
	case REALSXP:
		_IRanges_memcpy_to_i1i2(i1, i2,
			(char *) REAL(out), LENGTH(out),
			(char *) REAL(in), LENGTH(in), sizeof(double));
		break;
	case CPLXSXP:
		_IRanges_memcpy_to_i1i2(i1, i2,
			(char *) COMPLEX(out), LENGTH(out),
			(char *) COMPLEX(in), LENGTH(in), sizeof(Rcomplex));
		break;
	default:
		error("IRanges internal error in _vector_copy_to_offset(): "
		      "%s type not supported", type2str(TYPEOF(out)));
	}
	return;
}

/*
 * RECYCLING: Cyclic reading from 'in'.
 * INTERFACE: "subset".
 * In addition, "raw" vectors support fast on-the-fly translation via the
 * 'lkup' table.
 */
void _vector_copy_to_subset(SEXP out, SEXP in, SEXP subset, SEXP lkup)
{
	switch (TYPEOF(out)) {
	case RAWSXP:
		if (lkup == R_NilValue)
			_IRanges_memcpy_to_subset(
				INTEGER(subset), LENGTH(subset),
				(char *) RAW(out), LENGTH(out),
				(char *) RAW(in), LENGTH(in), sizeof(Rbyte));
		else
			_IRanges_charcpy_to_subset_with_lkup(
				INTEGER(subset), LENGTH(subset),
				(char *) RAW(out), LENGTH(out),
				(char *) RAW(in), LENGTH(in), 
				INTEGER(lkup), LENGTH(lkup));
		break;
	case LGLSXP:
	case INTSXP:
		_IRanges_memcpy_to_subset(INTEGER(subset), LENGTH(subset),
			(char *) INTEGER(out), LENGTH(out),
			(char *) INTEGER(in), LENGTH(in), sizeof(int));
		break;
	case REALSXP:
		_IRanges_memcpy_to_subset(INTEGER(subset), LENGTH(subset),
			(char *) REAL(out), LENGTH(out),
			(char *) REAL(in), LENGTH(in), sizeof(double));
		break;
	case CPLXSXP:
		_IRanges_memcpy_to_subset(INTEGER(subset), LENGTH(subset),
			(char *) COMPLEX(out), LENGTH(out),
			(char *) COMPLEX(in), LENGTH(in), sizeof(Rcomplex));
		break;
	default:
		error("IRanges internal error in _vector_copy_to_subset(): "
		      "%s type not supported", type2str(TYPEOF(out)));
	}
	return;
}
