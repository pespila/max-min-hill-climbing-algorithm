// #include <Rcpp.h>
#include <omp.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <vector>
#include <R_ext/Arith.h>
#include <R_ext/Utils.h>
#include <tr1/unordered_map>
using namespace std;

int* MyTest() {
	std::tr1::unordered_map<int, int> d;
	d[1] = 2;
	int *x;
	x = (int*)R_alloc(5, sizeof(int));
	printf("%d\n", x[1]);
	return x;
}