#include <omp.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/Arith.h>
#include <R_ext/Utils.h>

SEXP MyTest() {
	// std::tr1::unordered_map<int, int> d;
	// d[1] = 2;
	SEXP x;
	PROTECT(x = allocVector(INTSXP,5));
	// int *x;
	// x = (int*)R_alloc(5, sizeof(int));
	// printf("%d\n", x[1]);
	UNPROTECT(1);
	return x;
}