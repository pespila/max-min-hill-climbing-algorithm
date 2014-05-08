#include <Rcpp.h>
#include <omp.h>
#include <Rinternals.h>
#include <R.h>
using namespace Rcpp;

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

bool allC(SEXP a, SEXP b, int del = -1) {
	bool out = TRUE;
	NumericVector x(a), y(b);
	if(x.size()!=y.size()) out=FALSE;
	else {
		for (int i = 0; i < x.size(); i++) {
			if(i == del) continue;
			if(del == -2 && (i == 0 || i == 1)) continue;
			if(x[i]!=y[i]) {out=FALSE;break;}
		}
	}
	return out;
}

// [[Rcpp::export]]
NumericVector Cardinality(SEXP x) {
	NumericMatrix A(x);
	NumericVector out;
	int m = A.nrow(), n = A.ncol();

	for (int i = 0; i < n; i++)
	{
		NumericVector count;
		for (int j = 0; j < m; j++)
		{
			if(i == m)
				count.push_back(A(i,j));
			for (int k = 0; k < count.size(); k++)
			{
				if(A(i,j) != count[k])
					count.push_back(A(i,j));
			}
		}
		out.push_back(count.size());
	}
	return out;
}

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
double Statistic(SEXP x, SEXP y) {
	int n, m;
	NumericVector count(4,0.0);
	NumericVector sum(1,0.0), pvalue(1,0.0);
	NumericVector df;
	
	NumericMatrix A(x);
	NumericMatrix B(y);
	n = A.nrow();
	m = B.nrow();
	NumericVector u;
	NumericVector v;

	#pragma omp parallel
	{
		#pragma omp for 
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

	df = Cardinality(A);
	
	int DF = (df[0] - 1) * (df[1] - 1);

	for (int i = 2; i < df.size(); i++)
	{
		DF *= df[i];
	}

	pvalue = pchisq(sum, DF, FALSE);
	return pvalue[0];
}

// [[Rcpp::export]]
int S(SEXP a,SEXP X,int type=1) {
	int out = 0;
	if(type) {
		NumericVector x(X);
		NumericMatrix A(a);
		NumericVector tmp(A.nrow());

		int n=A.nrow(),m=A.ncol();
		#pragma omp parallel
		{
			#pragma omp for
			for (int i = 0; i < n; i++)
			{
				tmp=A.row(i);
				if(allC(x,tmp))
					out++;
			}
		}
	}
	else {
		CharacterVector x(X);
		CharacterMatrix A(a);
		CharacterVector tmp(A.nrow());

		int n=A.nrow(),m=A.ncol();
		#pragma omp parallel
		{
			#pragma omp for
			for (int i = 0; i < n; i++)
			{
				tmp=A.row(i);
				if(allC(x,tmp))
					out++;
			}
		}
	}
	return out;
}

// // [[Rcpp::export]]
// NumericVector Run(SEXP x,SEXP y) {
// 	CharacterMatrix A(x);
// 	CharacterMatrix B(y);
// 	CharacterVector tmpOrig;
// 	CharacterVector tmpUniqu;
// 	NumericVector out(B.nrow());

// 	#pragma omp parallel
// 	{
// 		#pragma omp for
// 		for (int i = 0; i < B.nrow(); i++)
// 		{
// 			for (int j = 0; j < A.nrow(); j++)
// 			{
// 				tmpOrig = A.row(j);
// 				tmpUniqu = B.row(i);
// 				if (allCp(tmpOrig, tmpUniqu))
// 					out[i]++;
// 			}
// 		}
// 	}
// 	return out;
// }

// [[Rcpp::export]]
int Df(NumericVector DfSet) {
	int df;
	df=(DfSet[0]-1)*(DfSet[1]-1);
	#pragma omp parallel
	{
		#pragma omp for 
		for (int i = 2; i < DfSet.size(); i++)
		{
			df*=DfSet[i];
		}
	}
	return df;
}