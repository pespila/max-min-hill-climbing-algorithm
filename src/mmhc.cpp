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