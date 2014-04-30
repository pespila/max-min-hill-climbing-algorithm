#include <Rcpp.h>
using namespace Rcpp;

class Main
{
public:
	Main();
	~Main();

	void Shout();
};

Main::Main() {}

Main::~Main() {}

void Main::Shout() {
	Rprintf("Yippie\n");
}

// [[Rcpp::export]]
NumericVector convolveCpp(NumericVector a, NumericVector b) {
	int na = a.size(), nb = b.size();
	int nab = na + nb - 1;
	NumericVector xab(nab);
	for (int i = 0; i < na; i++)
		for (int j = 0; j < nb; j++)
			xab[i + j] += a[i] * b[j];
	return xab;
}

// [[Rcpp::export]]
bool equal(NumericVector x, NumericVector y) {
	int n=0;
	while(x[n]==y[n]) {
		if(n==x.size()) return TRUE;
		n++;
	}
	return FALSE;
}

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

// [[Rcpp::export]]
NumericVector MVM(NumericMatrix A, NumericVector& x) {
	NumericVector Ax(x.size());
	double sum=0.0;
	for (int i = 0; i < x.size(); i++) {
		for (int j = 0; j < x.size(); j++) {
			Ax[i]+=A(i,j)*x[j];
		}
	}
	Main B;
	B.Shout();
	return Ax;
	// x=Ax;


	// return Ax;
}

// // [[Rcpp:export]]
// int Sabc(NumericVector x, NumericMatrix A) {
// 	int i=0,j=0,n=A.nrow(),m=A.ncol(),S_abc=0,S_abc_tmp=0;
// 	for (int i; i < n; i++)
// 	{
// 		S_abc_tmp=0;
// 		for (int j; j < m; j++)
// 		{
// 			if(x[j]==A(i,j)) S_abc_tmp++;
// 		}
// 		if(S_abc_tmp==m) S_abc++;
// 	}
// 	return S_abc;
// }

// [[Rcpp::export]]
double innerProduct(NumericVector x) {
	// double x = 0;
	return mean(x);
	// return rnorm(1, 1 / (x + 1), 1 / sqrt(2 * (x + 1)));
	// NumericVector Ax(x.size());
	// double sum=0.0;
	// for (int i = 0; i < x.size(); i++) {
	// 	for (int j = 0; j < x.size(); j++) {
	// 		for (int k = 0; k < x.size(); k++) {
	// 			sum+=A(j,k)*x[k];
	// 		}
	// 	}
	// 	Ax[i]=exp(sum);
	// 	sum=0.0;
	// }
	// x=Ax;


	// return Ax;
}