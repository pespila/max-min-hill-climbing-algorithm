#include <Rcpp.h>
#include <tr1/unordered_map>

using namespace std;
using namespace std::tr1;
using namespace Rcpp;

// [[Rcpp::export]]
SEXP test(SEXP x) {
	DataFrame A(x);
	NumericMatrix B(A.nrows(), A.length());
	NumericVector u;
	for (int i = 0; i < A.length(); i++)
	{
		u = A[i];
		B(_, i) = u;
	}
	return B;
}