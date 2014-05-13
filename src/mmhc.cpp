#include <Rcpp.h>
#include <omp.h>
#include <Rinternals.h>
#include <R.h>
#include <vector>
using namespace Rcpp;
using namespace std;

// bool allC(SEXP a, SEXP b, int del = -1) {
// 	bool out = TRUE;
// 	switch (TYPEOF(a)) {
// 		case INTSXP:
// 		{
// 			IntegerVector x(a), y(b);
// 			if(x.size()!=y.size()) out=FALSE;
// 			else {
// 				for (int i = 0; i < x.size(); i++) {
// 					if(i == del) continue;
// 					if(del == -2 && (i == 0 || i == 1)) continue;
// 					if(x[i]!=y[i]) {out=FALSE;break;}
// 				}
// 			}
// 			break;
// 		}
// 		case REALSXP:
// 		{
// 			NumericVector x(a), y(b);
// 			if(x.size()!=y.size()) out=FALSE;
// 			else {
// 				for (int i = 0; i < x.size(); i++) {
// 					if(i == del) continue;
// 					if(del == -2 && (i == 0 || i == 1)) continue;
// 					if(x[i]!=y[i]) {out=FALSE;break;}
// 				}
// 			}
// 			break;
// 		}
// 		default:
// 		{
// 			CharacterVector x(a), y(b);
// 			if(x.size()!=y.size()) out=FALSE;
// 			else {
// 				for (int i = 0; i < x.size(); i++) {
// 					if(i == del) continue;
// 					if(del == -2 && (i == 0 || i == 1)) continue;
// 					if(x[i]!=y[i]) {out=FALSE;break;}
// 				}
// 			}
// 			break;
// 		}
// 	}
// 	return out;
// }

// // [[Rcpp::export]]
// NumericVector Cardinality(SEXP x) {
// 	NumericMatrix A(x);
// 	NumericVector out;
// 	int m = A.nrow(), n = A.ncol();
// 	vector<double> count;

// 	#pragma omp parallel
// 	{
// 		#pragma omp for collapse(2)
// 		for (int i = 0; i < n; i++)
// 		{
// 			for (int j = 0; j < m; j++)
// 			{
// 				count.push_back(A(j,i));
// 				for (int k = 0; k < count.size(); k++)
// 				{
// 					if(A(j,i) != count[k])
// 						count.push_back(A(j,i));
// 				}
// 			}
// 			out.push_back(count.size());
// 			count.clear();
// 		}
// 	}
// 	return out;
// }

// // [[Rcpp::export]]
// double Statistic(SEXP x, SEXP y, SEXP z) {
// 	int n, m;
// 	NumericVector count(4,0.0);
// 	NumericVector sum(1,0.0), pvalue(1,0.0);
// 	NumericVector df(z);

// 	switch(TYPEOF(x)) {
// 		case INTSXP:
// 		{
// 			IntegerMatrix A(x);
// 			IntegerMatrix B(y);
// 			n = A.nrow();
// 			m = B.nrow();
// 			IntegerVector u;
// 			IntegerVector v;

// 			#pragma omp parallel
// 			{
// 				#pragma omp for 
// 				for (int i = 0; i < m; i++)
// 				{
// 					count[0] = count[1] = count[2] = count[3] = 0;
// 					u = B.row(i);
// 					for (int j = 0; j < n; j++)
// 					{
// 						v = A.row(j);
// 						if(allC(v,u))
// 							count[0]++;
// 						if(allC(v,u,1))
// 							count[1]++;
// 						if(allC(v,u,0))
// 							count[2]++;
// 						if(v.size() == n) 
// 							count[3] = n;
// 						else {
// 							if(allC(v,u,-2))
// 								count[3]++;
// 						}
// 					}
// 					sum[0] += 2 * count[0] * log( ( count[0] * count[3] ) / ( count[1] * count[2] ) );
// 				}
// 			}
// 			break;
// 		}
// 		case REALSXP:
// 		{
// 			NumericMatrix A(x);
// 			NumericMatrix B(y);
// 			n = A.nrow();
// 			m = B.nrow();
// 			NumericVector u;
// 			NumericVector v;

// 			#pragma omp parallel
// 			{
// 				#pragma omp for 
// 				for (int i = 0; i < m; i++)
// 				{
// 					count[0] = count[1] = count[2] = count[3] = 0;
// 					u = B.row(i);
// 					for (int j = 0; j < n; j++)
// 					{
// 						v = A.row(j);
// 						if(allC(v,u))
// 							count[0]++;
// 						if(allC(v,u,1))
// 							count[1]++;
// 						if(allC(v,u,0))
// 							count[2]++;
// 						if(allC(v,u,-2))
// 							count[3]++;
// 					}
// 					sum[0] += 2 * count[0] * log( ( count[0] * count[3] ) / ( count[1] * count[2] ) );
// 				}
// 			}
// 			break;
// 		}
// 		default:
// 		{
// 			CharacterMatrix A(x);
// 			CharacterMatrix B(y);
// 			n = A.nrow();
// 			m = B.nrow();
// 			CharacterVector u;
// 			CharacterVector v;

// 			#pragma omp parallel
// 			{
// 				#pragma omp for 
// 				for (int i = 0; i < m; i++)
// 				{
// 					count[0] = count[1] = count[2] = count[3] = 0;
// 					u = B.row(i);
// 					for (int j = 0; j < n; j++)
// 					{
// 						v = A.row(j);
// 						if(allC(v,u))
// 							count[0]++;
// 						if(allC(v,u,1))
// 							count[1]++;
// 						if(allC(v,u,0))
// 							count[2]++;
// 						if(allC(v,u,-2))
// 							count[3]++;
// 					}
// 					sum[0] += 2 * count[0] * log( ( count[0] * count[3] ) / ( count[1] * count[2] ) );
// 				}
// 			}
// 			break;
// 		}
// 	}

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

// // [[Rcpp::export]]
// SEXP Unique(SEXP x) {
// 	NumericMatrix A(x);
// 	NumericVector u(A.ncol()), w(A.ncol());
// 	int cols = A.ncol();
// 	int** U = new int*[cols];
// 	for (int i = 0; i < rows; ++i)
// 	NumericMatrix U;
// 	int n = 0;
// 	for (int i = 0; i < A.nrow(); i++)
// 	{
//     	U[i] = new int[cols];
// 		v = A(i, _);
// 		if (All(u, w)) {
// 			continue;
// 		} else {
// 			for (int j = 0; j < U.ncol(); j++)
// 			{
// 				U(n, j) = A(i, j);
// 			}
// 			// U(n, _) = u;
// 			n++;
// 		}
// 	}
// 	return A;
// }

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