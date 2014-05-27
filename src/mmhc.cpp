#include <RcppCommon.h>
#include <Rcpp.h>
#include <omp.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <vector>
#include <R_ext/Arith.h>
#include <R_ext/Utils.h>
// #include <unordered_map>
#include <tr1/unordered_map>
// using namespace Rcpp;
using namespace std;
using namespace std::tr1;

using namespace Rcpp;


class my
{
private:
	NumericVector x;
public:
	my(int n) {}
	~my();
	NumericVector Get() {
		return this->x;
	}
};

// namespace Rcpp {
// 	template <my> SEXP wrap( const my& ) ;
// }


/// create an external pointer to a Uniform object
RcppExport SEXP my__new(SEXP n) {
	// convert inputs to appropriate C++ types
	NumericVector x(as<int>(n),1.0);
	// int N = as<int>(n);
	// create a pointer to an Uniform object and wrap it
	// as an external pointer
	Rcpp::XPtr<my> ptr( new my( as<int>(n) ), true );
	// return the external pointer to the R side
	return ptr;
}

/// invoke the draw method
RcppExport SEXP my__Get( SEXP xp ) {
	// grab the object as a XPtr (smart pointer) to Uniform
	Rcpp::XPtr<my> ptr(xp);
	// convert the parameter to int
	// invoke the function
	NumericVector res = ptr->Get();
	// return the result to R
	return res;
}

// [[plugins(cpp11)]]

// 	int DF = (df[0] - 1) * (df[1] - 1);

// 	for (int i = 2; i < df.size(); i++)
// 	{
// 		DF *= df[i];
// 	}

// 	pvalue = pchisq(sum, DF, FALSE);
// 	return pvalue[0];
// }


// [[Rcpp::export]]
bool allC(SEXP a, SEXP b, int del = -1) {
	bool out = TRUE;
	double *x = REAL(a), *y = REAL(b);
	if(XLENGTH(a)!=XLENGTH(b)) out=FALSE;
	else {
		for (int i = 0; i < XLENGTH(a); i++) {
			if(i == del) continue;
			if(del == -2 && XLENGTH(b) == 2) {out=TRUE; break;}
			if(del == -2 && (i == 0 || i == 1)) continue;
			if(x[i]!=y[i]) {
				out=FALSE;
				break;
			}
		}
	}
	return out;
}

// [[Rcpp::export]]
bool All(SEXP a, SEXP b) {
	bool out = FALSE;
	NumericMatrix A(a);
	double *x = REAL(b);
	// NumericVector x(b);
	int n;

	for (int i = 0; i < A.nrow(); i++)
	{
		n = 0;
		for (int j = 0; j < A.ncol(); j++)
		{
			if (A(i, j) == x[j]) {
				n++;
			} else {
				break;
			}
		}
		if (n == A.ncol()) {
			out = TRUE;
			break;
		}
	}
	return out;
}

// [[Rcpp::export]]
int Df(SEXP x) {
	double *DfSet = REAL(x);
	int df = (DfSet[0]-1)*(DfSet[1]-1);
	for (int i = 2; i < XLENGTH(x); i++)
	{
		df*=DfSet[i];
	}
	return df;
}

// [[Rcpp::export]]
NumericVector Statistic(SEXP x, SEXP y, SEXP z) {
	int n, m;
	NumericVector count(4,0.0);
	NumericVector sum(1,0.0);
	NumericVector pvalue(1,0.0);
	NumericVector df(z);
	NumericVector out(2, 0.0);
	
	NumericMatrix A(x);
	NumericMatrix B(y);
	n = A.nrow();
	m = B.nrow();
	NumericVector u;
	NumericVector v;

	#pragma omp parallel
	{
		#pragma omp for collapse(2)
		for (int i = 0; i < m; i++)
		{
			count[0] = count[1] = count[2] = count[3] = 0;
			u = B.row(i);
			for (int j = 0; j < n; j++)
			{
				v = A.row(j);
				if(allC(v,u))
					count[0]++;
				if(allC(v,u,1))
					count[1]++;
				if(allC(v,u,0))
					count[2]++;
				if(allC(v,u,-2))
					count[3]++;
			}
			sum[0] += 2 * count[0] * log( ( count[0] * count[3] ) / ( count[1] * count[2] ) );
		}
	}

	// df = Cardinality(A);
	int DF = Df(df);
	pvalue = pchisq(sum, DF, FALSE);
	out[0] = pvalue[0];
	out[1] = sum[0];

	return out;
}

// // [[Rcpp::export]]
// double Norm(SEXP x, SEXP y, SEXP a) {
// 	IntegerVector X(x);
// 	IntegerVector Y(y);
// 	IntegerMatrix A(a);
// 	int n = X.size(), i = 0, j = 0, I = A.ncol(), J = A.nrow();
// 	IntegerVector tmp(n,0);

// 	double out=0;

// 	for (i = 0; i < I; i++)
// 	{
// 		for (j = 0; j < J; j++)
// 		{
// 			tmp[j] += A(j,i) * X[i];
// 		}
// 	}

// 	for (i = 0; i < n; i++)
// 	{
// 		out += tmp[i] + Y[i];
// 	}
// 	out = sqrt(out);
// 	return out;
// }

// // [[Rcpp::export]]
// double NoRm(SEXP x, SEXP y, SEXP a, SEXP k, SEXP l) {
// 	int *X = INTEGER(x);
// 	int *Y = INTEGER(y);
// 	int *A = INTEGER(a);
// 	int I = INTEGER(k)[0];
// 	int J = INTEGER(l)[0];
// 	int n = XLENGTH(x), i = 0, j = 0;
// 	int *tmp = (int*)R_alloc(n, sizeof(int));
	
// 	double out=0;

// 	for (i = 0; i < n; i++)
// 	{
// 		tmp[i] = 0;
// 	}

// 	for (i = 0; i < I; i++)
// 	{
// 		for (j = 0; j < J; j++)
// 		{
// 			tmp[j] += A[j + I * i]*X[i];
// 		}
// 	}

// 	for (i = 0; i < n; i++)
// 	{
// 		out += tmp[i] + Y[i];
// 	}

// 	out = sqrt(out);
// 	return out;
// }

// unordered_map<NumericVector, int> Map(SEXP m) {
// 	unordered_map<int*, int> mapping;
// 	NumericMatrix A(m);
// 	int row = A.nrow(), col = A.ncol();

// 	double *tmp = (double*)R_alloc(col, sizeof(double));

// 	for (int i = 0; i < row; i++)
// 	{
// 		tmp = A(i, _);
// 	}

// 	// for (int i = 0; i < col; i++)
// 	// {
// 	// 	for (int j = 0; j < col; j++)
// 	// 	{
// 	// 		tmp[j] = matrix[i + row * j];
// 	// 	}
// 	// 	mapping[tmp] = 1;
// 	// }
// 	return mapping;
// }

// unordered_map<NumericVector, int> Map(SEXP m, int col) {
// 	unordered_map<NumericVector, int> mapping;
// 	double *matrix = REAL(m);
// 	int row = XLENGTH(m)/col;

// 	NumericVector tmp(col);

// 	for (int i = 0; i < row; i++)
// 	{
// 		for (int j = 0; j < col; j++)
// 		{
// 			tmp[j] = matrix[i + row * j];
// 		}
// 		mapping[tmp] = 1;
// 	}
// 	return mapping;
// }

bool In(int x, NumericVector& y) {

	bool out = FALSE;
	for (int i = 0; i < y.size(); i++)
	{
		if (y[i] == x) {
			out = TRUE;
			break;
		}
	}

	return out;

}

// [[Rcpp::export]]
SEXP Calc(SEXP m) {
	NumericMatrix A(m);

	for (int i = 0; i < A.ncol(); i++)
	{
		for (int j = 0; j < A.nrow(); j++)
		{
			A(j,i) = sqrt(A(j,i) / sqrt(3));
		}
	}

	return A;
}