#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void Gsquare(NumericVector a) {
	for (int i = 0; i < a.size(); i++)
	{
		Rprintf("%f ",a[i]);
	}
	Rprintf("\n");
}

// [[Rcpp::export]]
int Sabc(NumericVector x, NumericMatrix A) {
	int n=A.nrow(),m=A.ncol(),S_abc=0,S_abc_tmp=0;
	for (int i = 0; i < n; i++)
	{
		S_abc_tmp=0;
		for (int j = 0; j < m; j++)
		{
			if(x[j]==A(i,j)) S_abc_tmp++;
		}
		if(S_abc_tmp==m) S_abc++;
	}
	return S_abc;
	// return x.size();
}