// #include <RcppCommon.h>
// #include <Rcpp.h>
#include <omp.h>
#include <R.h>
// #include <Rinternals.h>
#include <Rmath.h>
// #include <Rdefines.h>
#include <vector>
// #include <R_ext/Arith.h>
// #include <R_ext/Utils.h>
// #include <unordered_map>
// #include <boost/lexical_cast.hpp>
// #include <boost>
#include <tr1/unordered_map>
#include <RcppArmadillo.h>
#include <time.h>
#include <string>
#include <sstream>
using namespace std;
using namespace std::tr1;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
string HashC(IntegerVector& array, int i, bool skip) {
	ostringstream oss("");
	for (i; i < array.size(); i++) {
		if (skip && i == 1)
			continue;
		oss << array[i];
	}
	
	return oss.str();
}

// [[Rcpp::export]]
IntegerMatrix TP(IntegerMatrix& A, int row, int col) {
	IntegerMatrix B(row, col);
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			B(i, j) = A(j, i);
	return B;
}

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
	IntegerMatrix A(mat);
	int hDim = A.ncol(), vDim = A.nrow();
	NumericVector sum(1, 0.0), pvalue(1, 0.0), out(2, 0.0);
	IntegerVector cardinality(hDim);

	for (int i = 0; i < A.ncol(); i++)
	{
		cardinality[i] = max(A(_, i));
	}

	if (hDim == 2) {
		int l = cardinality[1], m = cardinality[0];
		double *x = OneD(m), *y = OneD(l), **z = TwoD(m, l);

		for (int i = 0; i < vDim; i++)
		{
			x[A(i, 0) - 1]++;
			y[A(i, 1) - 1]++;
			z[A(i, 0) - 1][A(i, 1) - 1]++;
		}

		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < l; j++)
			{
				if (x[i] < 1 || y[j] < 1 || z[i][j] < 1)
					continue;
				sum[0] += 2.0 * z[i][j] * log( (z[i][j] * vDim) / (x[i] * y[j]) );
				// cout << "2.0 * " << z[i][j] << "* log( ( " << z[i][j] << " * " << vDim << " ) / ( " << x[i] << " * " << y[j] << ") ) = " << 2.0 * z[i][j] * log( (z[i][j] * vDim) / (x[i] * y[j]) ) << endl;
			}
		}

		int DF = cardinality[0] * cardinality[1];
		for (int i = 2; i < cardinality.size(); i++)
		{
			DF *= (cardinality[i] + 1);
		}
		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;

	}  else if (hDim == 3) {
		int k = cardinality[2], l = cardinality[1], m = cardinality[0];
		double *v = OneD(k), **x = TwoD(m, k), **y = TwoD(l, k), ***z = ThreeD(m, l, k);

		for (int i = 0; i < vDim; i++)
		{
			v[A(i, 2) - 1]++;
			x[A(i, 0) - 1][A(i, 2) - 1]++;
			y[A(i, 1) - 1][A(i, 2) - 1]++;
			z[A(i, 0) - 1][A(i, 1) - 1][A(i, 2) - 1]++;
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

		int DF = cardinality[0] * cardinality[1];
		for (int i = 2; i < cardinality.size(); i++)
		{
			DF *= (cardinality[i] + 1);
		}
		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;

	} else if (hDim == 4) {
		int n = cardinality[3], k = cardinality[2], l = cardinality[1], m = cardinality[0];
		double **v = TwoD(k, n), ***x = ThreeD(m, k, n), ***y = ThreeD(l, k, n), ****z = FourD(m, l, k, n);

		for (int i = 0; i < vDim; i++)
		{
			v[A(i, 2) - 1][A(i, 3) - 1]++;
			x[A(i, 0) - 1][A(i, 2) - 1][A(i, 3) - 1]++;
			y[A(i, 1) - 1][A(i, 2) - 1][A(i, 3) - 1]++;
			z[A(i, 0) - 1][A(i, 1) - 1][A(i, 2) - 1][A(i, 3) - 1]++;
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

		int DF = cardinality[0] * cardinality[1];
		for (int i = 2; i < cardinality.size(); i++)
		{
			DF *= (cardinality[i] + 1);
		}
		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;

	} else if (hDim == 5) {
		int o = cardinality[4], n = cardinality[3], k = cardinality[2], l = cardinality[1], m = cardinality[0];
		double ***v = ThreeD(k, n, o), ****x = FourD(m, k, n, o), ****y = FourD(l, k, n, o), *****z = FiveD(m, l, k, n, o);

		for (int i = 0; i < vDim; i++)
		{
			v[A(i, 2) - 1][A(i, 3) - 1][A(i, 4) - 1]++;
			x[A(i, 0) - 1][A(i, 2) - 1][A(i, 3) - 1][A(i, 4) - 1]++;
			y[A(i, 1) - 1][A(i, 2) - 1][A(i, 3) - 1][A(i, 4) - 1]++;
			z[A(i, 0) - 1][A(i, 1) - 1][A(i, 2) - 1][A(i, 3) - 1][A(i, 4) - 1]++;
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

		int DF = cardinality[0] * cardinality[1];
		for (int i = 2; i < cardinality.size(); i++)
		{
			DF *= (cardinality[i] + 1);
		}
		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;
		
	} else {
		string acKey, bcKey, abcKey, cKey;
		IntegerMatrix B = TP(A, hDim, vDim);
		unordered_map<string, int> Count;
		unordered_map<int, string> ReMap;
		int teta = 0;
		IntegerVector tmp(hDim);

		for (int i = 0; i < vDim; i++)
		{
			tmp = B(_, i);
			abcKey = HashC(tmp, 0, FALSE);
			bcKey = HashC(tmp, 1, FALSE);
			cKey = HashC(tmp, 2, FALSE);
			acKey = HashC(tmp, 0, TRUE);
			if (Count[abcKey] == 0) {
				ReMap[teta] = abcKey;
				ReMap[teta+1] = acKey;
				ReMap[teta+2] = bcKey;
				ReMap[teta+3] = cKey;
				teta = teta + 4;
				Count[abcKey]++;
				Count[acKey]++;
				Count[bcKey]++;
				Count[cKey]++;
			}
		}

		for (int i = 0; i < ReMap.size(); i = i + 4)
		{
			sum[0] += 2.0 * Count[ReMap[i]] == 0 * log( (Count[ReMap[i]] == 0 * Count[ReMap[i+3]] == 0) / (Count[ReMap[i+1]] == 0 * Count[ReMap[i+2]] == 0) );
		}

		int DF = cardinality[0] * cardinality[1];
		for (int i = 2; i < cardinality.size(); i++)
		{
			DF *= (cardinality[i] + 1);
		}
		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;

	}
}

// [[Rcpp::export]]
SEXP MySvalue(SEXP mat) {
	IntegerMatrix A(mat);
	int hDim = A.ncol(), vDim = A.nrow();
	NumericVector sum(1, 0.0), pvalue(1, 0.0), out(2, 0.0);
	IntegerVector cardinality(hDim);

	for (int i = 0; i < A.ncol(); i++)
	{
		cardinality[i] = max(A(_, i));
	}

	if (hDim == 2) {
		int l = cardinality[1], m = cardinality[0];
		double *x = OneD(m), *y = OneD(l), **z = TwoD(m, l);

		for (int i = 0; i < vDim; i++)
		{
			x[A(i, 0) - 1]++;
			y[A(i, 1) - 1]++;
			z[A(i, 0) - 1][A(i, 1) - 1]++;
		}

		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < l; j++)
			{
				if (x[i] < 1 || y[j] < 1 || z[i][j] < 1)
					continue;
				sum[0] += 2.0 * z[i][j] * log( (z[i][j] * vDim) / (x[i] * y[j]) );
				// cout << "2.0 * " << z[i][j] << "* log( ( " << z[i][j] << " * " << vDim << " ) / ( " << x[i] << " * " << y[j] << ") ) = " << 2.0 * z[i][j] * log( (z[i][j] * vDim) / (x[i] * y[j]) ) << endl;
			}
		}

		int DF = cardinality[0] * cardinality[1];
		for (int i = 2; i < cardinality.size(); i++)
		{
			DF *= (cardinality[i] + 1);
		}
		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;

	}  else if (hDim == 3) {
		int k = cardinality[2], l = cardinality[1], m = cardinality[0];
		double *v = OneD(k), **x = TwoD(m, k), **y = TwoD(l, k), ***z = ThreeD(m, l, k);

		for (int i = 0; i < vDim; i++)
		{
			v[A(i, 2) - 1]++;
			x[A(i, 0) - 1][A(i, 2) - 1]++;
			y[A(i, 1) - 1][A(i, 2) - 1]++;
			z[A(i, 0) - 1][A(i, 1) - 1][A(i, 2) - 1]++;
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

		int DF = cardinality[0] * cardinality[1];
		for (int i = 2; i < cardinality.size(); i++)
		{
			DF *= (cardinality[i] + 1);
		}
		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;

	} else if (hDim == 4) {
		int n = cardinality[3], k = cardinality[2], l = cardinality[1], m = cardinality[0];
		double **v = TwoD(k, n), ***x = ThreeD(m, k, n), ***y = ThreeD(l, k, n), ****z = FourD(m, l, k, n);

		for (int i = 0; i < vDim; i++)
		{
			v[A(i, 2) - 1][A(i, 3) - 1]++;
			x[A(i, 0) - 1][A(i, 2) - 1][A(i, 3) - 1]++;
			y[A(i, 1) - 1][A(i, 2) - 1][A(i, 3) - 1]++;
			z[A(i, 0) - 1][A(i, 1) - 1][A(i, 2) - 1][A(i, 3) - 1]++;
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

		int DF = cardinality[0] * cardinality[1];
		for (int i = 2; i < cardinality.size(); i++)
		{
			DF *= (cardinality[i] + 1);
		}
		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;

	} else if (hDim == 5) {
		int o = cardinality[4], n = cardinality[3], k = cardinality[2], l = cardinality[1], m = cardinality[0];
		double ***v = ThreeD(k, n, o), ****x = FourD(m, k, n, o), ****y = FourD(l, k, n, o), *****z = FiveD(m, l, k, n, o);

		for (int i = 0; i < vDim; i++)
		{
			v[A(i, 2) - 1][A(i, 3) - 1][A(i, 4) - 1]++;
			x[A(i, 0) - 1][A(i, 2) - 1][A(i, 3) - 1][A(i, 4) - 1]++;
			y[A(i, 1) - 1][A(i, 2) - 1][A(i, 3) - 1][A(i, 4) - 1]++;
			z[A(i, 0) - 1][A(i, 1) - 1][A(i, 2) - 1][A(i, 3) - 1][A(i, 4) - 1]++;
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

		int DF = cardinality[0] * cardinality[1];
		for (int i = 2; i < cardinality.size(); i++)
		{
			DF *= (cardinality[i] + 1);
		}
		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;
		
	} else {
		string acKey, bcKey, abcKey, cKey;
		IntegerMatrix B = TP(A, hDim, vDim);
		unordered_map<string, int> Count;
		unordered_map<int, string> ReMap;
		int teta = 0;
		IntegerVector tmp(hDim);

		for (int i = 0; i < vDim; i++)
		{
			tmp = B(_, i);
			abcKey = HashC(tmp, 0, FALSE);
			bcKey = HashC(tmp, 1, FALSE);
			cKey = HashC(tmp, 2, FALSE);
			acKey = HashC(tmp, 0, TRUE);
			if (Count[abcKey] == 0) {
				ReMap[teta] = abcKey;
				ReMap[teta+1] = acKey;
				ReMap[teta+2] = bcKey;
				ReMap[teta+3] = cKey;
				teta = teta + 4;
				Count[abcKey]++;
				Count[acKey]++;
				Count[bcKey]++;
				Count[cKey]++;
			}
		}

		for (int i = 0; i < ReMap.size(); i = i + 4)
		{
			sum[0] += 2.0 * Count[ReMap[i]] == 0 * log( (Count[ReMap[i]] == 0 * Count[ReMap[i+3]] == 0) / (Count[ReMap[i+1]] == 0 * Count[ReMap[i+2]] == 0) );
		}

		int DF = cardinality[0] * cardinality[1];
		for (int i = 2; i < cardinality.size(); i++)
		{
			DF *= (cardinality[i] + 1);
		}
		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;

	}
}

// [[Rcpp::export]]
SEXP TheValue(SEXP mat) {
	IntegerMatrix A(mat);
	int hDim = A.ncol(), vDim = A.nrow();
	IntegerMatrix B = TP(A, hDim, vDim);
	NumericVector sum(1, 0.0), pvalue(1, 0.0), out(2, 0.0);
	IntegerVector cardinality(A.ncol());

	for (int i = 0; i < A.ncol(); i++)
	{
		cardinality[i] = max(A(_, i));
	}

	if (hDim == 2) {
		int l = cardinality[1], m = cardinality[0];
		double *x = OneD(m), *y = OneD(l), **z = TwoD(m, l);

		for (int i = 0; i < vDim; i++)
		{
			x[A(i, 0) - 1]++;
			y[A(i, 1) - 1]++;
			z[A(i, 0) - 1][A(i, 1) - 1]++;
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

		int DF = m * l;
		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;

	}  else if (hDim == 3) {
		int k = cardinality[2], l = cardinality[1], m = cardinality[0];
		double *v = OneD(k), **x = TwoD(m, k), **y = TwoD(l, k), ***z = ThreeD(m, l, k);

		for (int i = 0; i < vDim; i++)
		{
			v[A(i, 2) - 1]++;
			x[A(i, 0) - 1][A(i, 2) - 1]++;
			y[A(i, 1) - 1][A(i, 2) - 1]++;
			z[A(i, 0) - 1][A(i, 1) - 1][A(i, 2) - 1]++;
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

		int DF = m * l * (k + 1);
		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;

	} else if (hDim > 3) {
		string acKey, bcKey, abcKey, cKey;
		unordered_map<string, int> Count;
		unordered_map<int, string> ReMap;
		int teta = 0;
		IntegerVector tmp(hDim);

		for (int i = 0; i < vDim; i++)
		{
			tmp = B(_, i);
			// tmp = A(i, _);
			abcKey = HashC(tmp, 0, FALSE);
			bcKey = HashC(tmp, 1, FALSE);
			cKey = HashC(tmp, 2, FALSE);
			acKey = HashC(tmp, 0, TRUE);
			if (Count[abcKey] == 0) {
				ReMap[teta] = abcKey;
				ReMap[teta+1] = acKey;
				ReMap[teta+2] = bcKey;
				ReMap[teta+3] = cKey;
				teta = teta + 4;
				Count[abcKey]++;
				Count[acKey]++;
				Count[bcKey]++;
				Count[cKey]++;
			}
		}

		for (int i = 0; i < ReMap.size(); i = i + 4)
		{
			sum[0] += 2.0 * Count[ReMap[i]] == 0 * log( (Count[ReMap[i]] == 0 * Count[ReMap[i+3]] == 0) / (Count[ReMap[i+1]] == 0 * Count[ReMap[i+2]] == 0) );
		}

		int DF = cardinality[0] * cardinality[1];
		for (int i = 2; i < cardinality.size(); i++)
		{
			DF *= (cardinality[i] + 1);
		}
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

void PrintLine(NumericVector A) {
	for (int i = 0; i < A.size(); i++)
	{
		cout << A[i] << " ";
	}
}

string int_array_to_string(NumericVector array) { //  int int_array[], int size_of_array){
	ostringstream oss("");
	for (int temp = 0; temp < array.size(); temp++)
		oss << array[temp];
	
	return oss.str();
}

// [[Rcpp::export]]
NumericMatrix Which(SEXP x) {
	srand(time(NULL));
	NumericMatrix A(x);
	NumericVector u;
	string key;
	unordered_map<string, int> Map;

	unordered_map<int, string> Uni;
	int n = Map.size(), k = 0;

	for (int i = 0; i < A.nrow(); i++)
	{
		u = A(i, _);
		key = int_array_to_string(u);
		Map[key] = rand()%10;
		if (Map.size() != n) {
			Uni[k] = key;
			n++;
			k++;
		}
	}

	NumericMatrix B(Map.size(), A.ncol());

	for (int i = 0; i < Map.size(); i++)
	{
		for (int j = 0; j < A.ncol(); j++)
		{
			B(i, j) = Uni[i][j] - '0';
		}
	}

	return B;
}

// [[Rcpp::export]]
IntegerMatrix Unique(SEXP x) {
	srand(time(NULL));
	IntegerMatrix A(x);
	NumericVector u;
	string key;
	unordered_map<string, int> Map;

	unordered_map<int, string> Uni;
	int n = Map.size(), k = 0;

	for (int i = 0; i < A.nrow(); i++)
	{
		u = A(i, _);
		key = int_array_to_string(u);
		Map[key] = rand()%10;
		if (Map.size() != n) {
			Uni[k] = key;
			n++;
			k++;
		}
	}

	IntegerMatrix B(Map.size(), A.ncol());

	for (int i = 0; i < Map.size(); i++)
	{
		for (int j = 0; j < A.ncol(); j++)
		{
			B(i, j) = Uni[i][j] - '0';
		}
	}

	return B;
}

// // [[Rcpp::export]]
// int Which(SEXP x, SEXP y) {
// 	int out = 0;
// 	NumericVector A(x);
// 	double *B = REAL(y);

// 	for (int i = 0; i < A.size(); i++)
// 	{
// 		if (A[i] == *B)
// 			out = i;
// 	}

// 	return out;
// }

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
NumericVector Statistics(SEXP x) {
	int n, m;
	NumericVector count(4,0.0);
	NumericVector sum(1,0.0);
	NumericVector pvalue(1,0.0);
	// NumericVector df(z);
	NumericVector out(2, 0.0);
	
	IntegerMatrix C(x);
	IntegerMatrix D = Unique(C);
	IntegerMatrix A = TP(C, C.ncol(), C.nrow());
	IntegerMatrix B = TP(D, D.ncol(), D.nrow());
	n = A.ncol();
	m = B.ncol();
	NumericVector u;
	NumericVector v;

	IntegerVector cardinality(A.nrow());

	for (int i = 0; i < A.nrow(); i++)
	{
		cardinality[i] = max(A(_, i));
	}

	for (int i = 0; i < m; i++)
	{
		count[0] = count[1] = count[2] = count[3] = 0;
		u = B(_, i);
		for (int j = 0; j < n; j++)
		{
			v = A(_, j);
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
		// cout << "2.0 * " << count[0] << " * log( " << count[0] << " * " << count[3] << " / " << count[1] << " * " << count[2] << " ) = " << 2.0 * count[0] * log( (count[0] * count[3]) / (count[1] * count[2]) ) << endl;
	}

	int DF = cardinality[0] * cardinality[1];
	for (int i = 2; i < cardinality.size(); i++)
	{
		DF *= (cardinality[i] + 1);
	}
	// int DF = Df(df);
	pvalue = pchisq(sum, DF, FALSE);
	out[0] = pvalue[0];
	out[1] = sum[0];

	return out;
}