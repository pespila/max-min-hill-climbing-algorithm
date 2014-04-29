#include <Rcpp.h>
#include <omp.h>
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

// [[Rcpp::export]]
void Gsquare(NumericVector x) {
	if(2) {
		Rprintf("Hi\n");
	}
}

// double Gsquare(NumericVector x, NumericMatrix A) {
// 	return 0.0;
// }

// NumericVector Sabc(NumericVector x, NumericMatrix A) {
// 	int n=A.nrow(),m=A.ncol(),S_abc=0,S_abc_tmp=0;
// 	NumericVector S;
// 	S.push_back(0);
// 	for (int i = 0; i < n; i++)
// 	{
// 		S_abc_tmp=0;
// 		for (int j = 0; j < m; j++)
// 		{
// 			if(x[j]==A(i,j)) S_abc_tmp++;
// 		}
// 		if(S_abc_tmp==m) {
// 			S[0]++;
// 			S.push_back(i);
// 		}
// 	}
// 	return S;
// 	// return x.size();
// }