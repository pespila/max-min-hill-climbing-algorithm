#include <Rcpp.h>
#include <omp.h>
#include <Rinternals.h>
#include <R.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Sabc(NumericMatrix A,NumericVector x) {
	int n=A.nrow(),m=A.ncol(),S_abc=0,S_abc_tmp=0;
	NumericVector out;
	out.push_back(0);
	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < n; i++)
		{
			S_abc_tmp=0;
			for (int j = 0; j < m; j++)
			{
				if(x[j]==A(i,j)) S_abc_tmp++;
			}
			if(S_abc_tmp==m) {
				out[0]++;
				out.push_back((-1)*(i+1));
			}
		}
	}
	return out;
}

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
DataFrame DaFr(SEXP data) {
	DataFrame frame(data);
	NumericVector out;
	out.push_back(0);
	return frame;
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
bool Proof(NumericVector x, NumericVector y) {
	// NumericVector out(x);
	// if(x==y) out=1;
	// else out=0;
	return is_true(any(x==y));
}

// [[Rcpp::export]]
int Df(NumericVector DfSet,int target,int selected,NumericVector subset) {
	int df;
	df=(DfSet[target-1]-1)*(DfSet[selected-1]-1);
	if(subset.size()==1) df*=DfSet[subset[0]-1];
	else {
		#pragma omp parallel
		{
			#pragma omp for 
			for (int i = 0; i < subset.size(); i++)
			{
				df*=DfSet[subset[i]-1];
			}
		}
	}
	return df;
}