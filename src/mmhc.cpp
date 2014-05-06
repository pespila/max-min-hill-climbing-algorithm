#include <Rcpp.h>
#include <omp.h>
#include <Rinternals.h>
#include <R.h>
using namespace Rcpp;

bool allCpp(NumericVector& x,NumericVector& y) {
	bool out=TRUE;
	if(x.size()!=y.size()) out=FALSE;
	else {
		for (int i = 0; i < x.size(); i++)
			if(x[i]!=y[i]) {out=FALSE;break;}
	}
	return out;
}

bool allCp(CharacterVector& x,CharacterVector& y) {
	bool out=TRUE;
	if(x.size()!=y.size()) out=FALSE;
	else {
		for (int i = 0; i < x.size(); i++)
			if(x[i]!=y[i]) {out=FALSE;break;}
	}
	return out;
}

// [[Rcpp::export]]
bool allC(SEXP s,SEXP t) {
	bool out=TRUE;
	int *x,*y;
	if(TYPEOF(s)==INTSXP) {x=INTEGER(s);y=INTEGER(t);}
	// if(TYPEOF(s)==REALSXP) {x=REAL(s);y=REAL(t);}
	// if(TYPEOF(s)==STRSXP) {x=CHAR(s);y=CHAR(t);}
	if(XLENGTH(s)!=XLENGTH(t)) out=FALSE;
	else {
		for (int i = 0; i < XLENGTH(s); i++)
			if(x[i]!=y[i]) {out=FALSE;break;}
	}
	return out;
}

// [[Rcpp::export]]
NumericVector Test(SEXP a,SEXP X,int type=1) {
	NumericVector out;
	out.push_back(0);
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
				if(allCpp(x,tmp)) {out[0]++; out.push_back((-1)*(i+1));}
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
				if(allCp(x,tmp)) {out[0]++; out.push_back((-1)*(i+1));}
			}
		}
	}
	return out;
}

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