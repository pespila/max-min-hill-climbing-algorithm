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
using namespace std;
using namespace std::tr1;
using namespace Rcpp;

// [[Rcpp::export]]
bool Id(SEXP a, SEXP b, int del = -1) {
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
SEXP UpdateCPC(SEXP x, int selected = 0) {
	NumericVector tmp;
	List cpc(x);
	if (cpc.size() == 0) {
		cpc.push_back(1);
		cpc.push_back(R_NilValue);
	} else if (cpc.size() == 2) {
		cpc[0] = 3;
		cpc.push_back(selected);
	} else {
		int cpcLength = cpc.size();
		cpc[0] = cpcLength + 1;
		cpc.push_back(selected);
		for (int i = 2; i < cpcLength; i++)
		{
			tmp = as<NumericVector>(cpc[i]);
			tmp.push_back(selected);
			cpc.push_back(tmp);
		}
	}

	return cpc;
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
NumericVector Statistics(SEXP x, SEXP y, SEXP z) {
	int n, m;
	NumericVector count(4, 0.0);
	NumericVector sum(1, 0.0);
	NumericVector pvalue(1, 0.0);
	NumericVector df(z);
	NumericVector out(2, 0.0);
	
	NumericMatrix A(x);
	NumericMatrix B(y);
	n = A.ncol();
	m = B.ncol();
	NumericVector u;
	NumericVector v;

	for (int i = 0; i < m; i++)
	{
		count[0] = count[1] = count[2] = count[3] = 0;
		u = B(_,i);
		for (int j = 0; j < n; j++)
		{
			v = A(_, j);
			if(Id(v, u))
				count[0]++;
			if(Id(v, u, 1))
				count[1]++;
			if(Id(v, u, 0))
				count[2]++;
			if(Id(v, u, -2))
				count[3]++;
		}
		sum[0] += 2 * count[0] * log( ( count[0] * count[3] ) / ( count[1] * count[2] ) );
	}

	// df = Cardinality(A);
	int DF = Df(df);
	pvalue = pchisq(sum, DF, FALSE);
	out[0] = pvalue[0];
	out[1] = sum[0];

	return out;
}