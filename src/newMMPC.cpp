// #include <RcppCommon.h>
// #include <Rcpp.h>
#include <omp.h>
#include <R.h>
// #include <Rinternals.h>
#include <Rmath.h>
// #include <Rdefines.h>
// #include <vector>
// #include <R_ext/Arith.h>
// #include <R_ext/Utils.h>
// #include <unordered_map>
// #include <tr1/unordered_map>
#include <RcppArmadillo.h>
using namespace std;
using namespace std::tr1;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(openmp)]]

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
SEXP T(SEXP x, SEXP n, SEXP m) {
	double *A = REAL(x);
	int *N = INTEGER(n), *M = INTEGER(m);
	NumericMatrix B(*N, *M);
	for (int i = 0; i < *M; i++)
		for (int j = 0; j < *N; j++)
			B(i, j) = A[j + i * (*M)];
	return B;
}

// [[Rcpp::export]]
void Parallel() {
	omp_set_num_threads(8);
	#pragma omp parallel
	cout << "Hi" << endl;
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

double *OneD(double x) {
	double *matrix = (double*)R_alloc(x, sizeof(double));
	for (int i = 0; i < x; i++)
	{
		matrix[i] = 0;
	}
	return matrix;
}

double **TwoD(double x, double y) {
	double **matrix = (double**)R_alloc(x, sizeof(double*));

	for (int i = 0; i < x; i++)
	{
		matrix[i] = (double*)R_alloc(y, sizeof(double));
		for (int j = 0; j < y; j++)
		{
			matrix[i][j] = 0;
		}
	}

	return matrix;
}

double ***ThreeD(double x, double y, double z) {
	double ***matrix = (double***)R_alloc(x, sizeof(double*));

	for (int i = 0; i < x; i++)
	{
		matrix[i] = (double**)R_alloc(y, sizeof(double*));
		for (int j = 0; j < y; j++)
		{
			matrix[i][j] = (double*)R_alloc(z, sizeof(double));
			for (int k = 0; k < z; k++)
			{
				matrix[i][j][k] = 0;
			}
		}
	}

	return matrix;
}

double ****FourD(double x, double y, double z, double a) {
	double ****matrix = (double****)R_alloc(x, sizeof(double*));

	for (int i = 0; i < x; i++)
	{
		matrix[i] = (double***)R_alloc(y, sizeof(double*));
		for (int j = 0; j < y; j++)
		{
			matrix[i][j] = (double**)R_alloc(z, sizeof(double*));
			for (int k = 0; k < z; k++)
			{
				matrix[i][j][k] = (double*)R_alloc(a, sizeof(double));
				for (int l = 0; l < a; l++)
				{
					matrix[i][j][k][l] = 0;
				}
			}
		}
	}

      // memset(p[i][j], '\0', sizeof(int) * depth);

	return matrix;
}

double *****FiveD(double x, double y, double z, double a, double b) {
	double *****matrix = (double*****)R_alloc(x, sizeof(double*));

	for (int i = 0; i < x; i++)
	{
		matrix[i] = (double****)R_alloc(y, sizeof(double*));
		for (int j = 0; j < y; j++)
		{
			matrix[i][j] = (double***)R_alloc(z, sizeof(double*));
			for (int k = 0; k < z; k++)
			{
				matrix[i][j][k] = (double**)R_alloc(a, sizeof(double*));
				for (int l = 0; l < a; l++)
				{
					matrix[i][j][k][l] = (double*)R_alloc(b, sizeof(double));
					for (int m = 0; m < b; m++)
					{
						matrix[i][j][k][l][m] = 0;
					}
				}
			}
		}
	}

	return matrix;
}

// [[Rcpp::export]]
SEXP MySvalue(SEXP mat) {
	NumericMatrix A(mat);
	int hDim = A.ncol(), vDim = A.nrow();
	NumericVector sum(1, 0.0), pvalue(1, 0.0), out(2, 0.0);

	// omp_set_num_threads(2);

	if (hDim == 2) {
		int l = max(A(_, 1)) + 1, m = max(A(_, 0)) + 1;
		double *x = OneD(m), *y = OneD(l), **z = TwoD(m, l);

		// #pragma omp parallel for
		for (int i = 0; i < vDim; i++)
		{
			x[(int)A(i, 0)]++;
			y[(int)A(i, 1)]++;
			z[(int)A(i, 0)][(int)A(i, 1)]++;
		}

		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < l; j++)
			{
				if (x[i] < 1 || y[j] < 1 || z[i][j] < 1)
					continue;
				sum[0] += 2.0 * z[i][j] * log( (z[i][j] * vDim) / (x[i] * y[j]) );
			}
		}

		int DF = (m-1) * (l-1);
		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;

	} else if (hDim == 3) {
		int k = max(A(_, 2)) + 1, l = max(A(_, 1)) + 1, m = max(A(_, 0)) + 1;
		double *v = OneD(k), **x = TwoD(m, k), **y = TwoD(l, k), ***z = ThreeD(m, l, k);

		// #pragma omp parallel for
		for (int i = 0; i < vDim; i++)
		{
			v[(int)A(i, 2)]++;
			x[(int)A(i, 0)][(int)A(i, 2)]++;
			y[(int)A(i, 1)][(int)A(i, 2)]++;
			z[(int)A(i, 0)][(int)A(i, 1)][(int)A(i, 2)]++;
		}

		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < l; j++)
			{
				for (int h = 0; h < k; h++)
				{
					if (x[i][h] < 1 || y[j][h] < 1 || z[i][j][h] < 1 || v[h] < 1)
						continue;
					sum[0] += 2.0 * z[i][j][h] * log( (z[i][j][h] * v[h]) / (x[i][h] * y[j][h]) );
				}
			}
		}

		int DF = (m-1) * (l-1) * k;
		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;

	} else if (hDim == 4) {
		int n = max(A(_, 3)) + 1, k = max(A(_, 2)) + 1, l = max(A(_, 1)) + 1, m = max(A(_, 0)) + 1;
		double **v = TwoD(k, n), ***x = ThreeD(m, k, n), ***y = ThreeD(l, k, n), ****z = FourD(m, l, k, n);

		// #pragma omp parallel for
		for (int i = 0; i < vDim; i++)
		{
			v[(int)A(i, 2)][(int)A(i, 3)]++;
			x[(int)A(i, 0)][(int)A(i, 2)][(int)A(i, 3)]++;
			y[(int)A(i, 1)][(int)A(i, 2)][(int)A(i, 3)]++;
			z[(int)A(i, 0)][(int)A(i, 1)][(int)A(i, 2)][(int)A(i, 3)]++;
		}

		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < l; j++)
			{
				for (int h = 0; h < k; h++)
				{
					for (int f = 0; f < n; f++)
					{
						if (x[i][h][f] < 1 || y[j][h][f] < 1 || z[i][j][h][f] < 1 || v[h][f] < 1)
							continue;
						sum[0] += 2.0 * z[i][j][h][f] * log( (z[i][j][h][f] * v[h][f]) / (x[i][h][f] * y[j][h][f]) );
					}
				}
			}
		}

		int DF = (m-1) * (l-1) * k * (n-1);
		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;

	} else if (hDim == 5) {
		int o = max(A(_, 4)) + 1, n = max(A(_, 3)) + 1, k = max(A(_, 2)) + 1, l = max(A(_, 1)) + 1, m = max(A(_, 0)) + 1;
		double ***v = ThreeD(k, n, o), ****x = FourD(m, k, n, o), ****y = FourD(l, k, n, o), *****z = FiveD(m, l, k, n, o);

		// #pragma omp parallel for
		for (int i = 0; i < vDim; i++)
		{
			v[(int)A(i, 2)][(int)A(i, 3)][(int)A(i, 4)]++;
			x[(int)A(i, 0)][(int)A(i, 2)][(int)A(i, 3)][(int)A(i, 4)]++;
			y[(int)A(i, 1)][(int)A(i, 2)][(int)A(i, 3)][(int)A(i, 4)]++;
			z[(int)A(i, 0)][(int)A(i, 1)][(int)A(i, 2)][(int)A(i, 3)][(int)A(i, 4)]++;
		}

		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < l; j++)
			{
				for (int h = 0; h < k; h++)
				{
					for (int f = 0; f < n; f++)
					{
						for (int e = 0; e < o; e++)
						{
							if (x[i][h][f][e] < 1 || y[j][h][f][e] < 1 || z[i][j][h][f][e] < 1 || v[h][f][e] < 1)
								continue;
							sum[0] += 2.0 * z[i][j][h][f][e] * log( (z[i][j][h][f][e] * v[h][f][e]) / (x[i][h][f][e] * y[j][h][f][e]) );
						}
					}
				}
			}
		}

		int DF = (m-1) * (l-1) * k * (n-1) * o;
		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;
		
	} else {
		out[0] = 1.0;
		out[1] = 1.0;

		return out;
	}
}

// [[Rcpp::export]]
int Df(SEXP x) {
	double *DfSet = REAL(x);
	int df = (DfSet[0]-1)*(DfSet[1]-1);
	// #pragma omp parallel for
	for (int i = 2; i < XLENGTH(x); i++)
	{
		df*=DfSet[i];
	}
	return df;
}

// [[Rcpp::export]]
NumericVector Statistics(SEXP x, SEXP y, SEXP z) {
	int n, m;
	NumericVector count(4,0.0);
	NumericVector sum(1,0.0);
	NumericVector pvalue(1,0.0);
	NumericVector df(z);
	NumericVector out(2, 0.0);
	
	NumericMatrix A(x);
	NumericMatrix B(y);
	n = A.ncol();
	m = B.ncol();
	NumericVector u;
	NumericVector v;

	#pragma omp parallel
	{
		#pragma omp for
		for (int i = 0; i < m; i++)
		{
			count[0] = count[1] = count[2] = count[3] = 0;
			u = B(_,i);
			for (int j = 0; j < n; j++)
			{
				v = A(_,j);
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

	int DF = Df(df);
	pvalue = pchisq(sum, DF, FALSE);
	out[0] = pvalue[0];
	out[1] = sum[0];

	return out;
}