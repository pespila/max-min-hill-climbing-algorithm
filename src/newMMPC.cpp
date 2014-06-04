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

// [[Rcpp::export]]
bool allC(SEXP a, SEXP b, int del = -1) {
	bool out = TRUE;
	double *x = REAL(a), *y = REAL(b);
	if(XLENGTH(a)!=XLENGTH(b)) out=FALSE;
	else {
		#pragma omp parallel
		{
			#pragma omp for
			{
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
		}
	}
	return out;
}

// [[Rcpp::export]]
SEXP UN(SEXP x) {
	int *X = INTEGER(x);
	typedef unordered_map<int,int> map_t;
	map_t Map;
	for (int i = 0; i < XLENGTH(x); i++)
	{
		Map.insert(map_t::value_type(X[i], 0));
		Map[X[i]]++;
	}
	return wrap(Map);
}

// double combine(double a, double b) {
//    	double times = 1;
//   //  	while (times <= b)
// 		// times *= 10;
//    	// return a*times + b;
//    	return a^b;
// } 

// // [[Rcpp::export]]
// SEXP UNB(SEXP x, SEXP n, SEXP m) {
// 	double *X = REAL(x);
// 	int *N = INTEGER(n), *M = INTEGER(m);
// 	double k;
// 	typedef unordered_map<int,int> map_t;
// 	map_t Map;
// 	a 
// 	for (int i = 0; i < *N; i++)
// 	{
// 		k = X[i];
// 		for (int j = 1; j < *M; j++)
// 		{
// 			k = combine(k, X[i + j * (*M)]);
// 		}
// 		Map.insert(map_t::value_type(k, 0));
// 		Map[k]++;
// 	}
// 	return wrap(Map);
// }

bool IsIn(mat X, mat R) {
	int row = X.n_rows, col = X.n_cols, test;
	bool equaled = FALSE;
	for (int i = 0; i < row; i++)
	{
		test = 0;
		for (int j = 0; j < col; j++)
		{
			if (X(i, j) == R(0, j))
				test++;
		}
		if (test == col) {
			equaled = TRUE;
			break;
		}
	}
	return equaled;
}

void insert(const mat& R, mat& A, int j) {
	for (int i = 0; i < R.n_cols; i++)
	{
		A(j, i) = R(0, i);
	}
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

// // [[Rcpp::export]]
// SEXP Test(SEXP x, SEXP n, SEXP m) {
// 	int *X = INTEGER(x), *N = INTEGER(n), *M = INTEGER(m);
// 	int K = (*N) * (*M);
// 	// int **A = (int**)R_alloc((*N), sizeof(int*));
// 	int *A = (int*)R_alloc(K, sizeof(int));
// 	cout << XLENGTH(A) << endl;
// 	for (int i = 0; i < K; i++)
// 	{
// 		// A[i] = X[i];
// 	}
// 	// for (int i = 0; i < *M; i++)
// 	// {
// 	// 	// A[i] = (int*)R_alloc((*M), sizeof(int));
// 	// 	for (int j = 0; j < *N; j++)
// 	// 	{
// 	// 		cout << X[j + i * (*M)] << endl;
// 	// 		// A[j + i * (*M)] = X[j + i * (*M)];
// 	// 	}
// 	// }
// 	return (SEXP)A;
// }

// // [[Rcpp::export]]
// mat Unique(SEXP x) {
// 	// if (A.n_cols >= A.n_rows) {
// 	// 	continue;
// 	// } else {
// 	// 	x = T(x);
// 	// }




// 	mat U(1, A.n_cols, fill::zeros);
// 	for (int i = 0; i < A.n_rows; i++)
// 	{
// 		if (!(IsIn(U, A.row(i)))) {
// 			insert(A.row(i), U, i);
// 			U.resize(U.n_rows+1, U.n_cols);
// 		}
// 	}
// 	return U;
// }

// // [[Rcpp::export]]
// mat MyTest(mat A) {
// 	return unique(A);
// }

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
		#pragma omp parallel
		{
			#pragma omp for
			{
				for (int i = 2; i < cpcLength; i++)
				{
					tmp = as<NumericVector>(cpc[i]);
					tmp.push_back(selected);
					cpc.push_back(tmp);
				}
			}
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

// // [[Rcpp::export]]
// SEXP Svalue(SEXP x) {
// 	NumericMatrix A(x);
// 	int m = max(A(_, 0)) + 1;

// 	if(A.ncol() == 1) {
// 		IntegerVector B(m, 0);

// 		for (int i = 0; i < A.nrow(); i++)
// 		{
// 			B[A[i]]++;
// 		}

// 		return B;
// 	} else if(A.ncol() == 2) {
// 		int l = max(A(_, 1)) + 1;
// 		NumericMatrix B(m, l);

// 		for (int i = 0; i < A.nrow(); i++)
// 		{
// 			B(A(i, 0), A(i, 1))++;
// 		}

// 		return B;
// 	} else if(A.ncol() == 3) {
// 		int l = max(A(_, 1)) + 1, k = max(A(_, 2)) + 1;
// 		double ***B = ThreeD(m, l, k);

// 		for (int i = 0; i < A.nrow(); i++)
// 		{
// 			B[(int)A(i, 0)][(int)A(i, 1)][(int)A(i, 2)]++;
// 		}

// 		return unique(A);
// 		// return B[0][0][0];

// 	}  else if(A.ncol() == 4) {
// 		int l = max(A(_, 1)) + 1, k = max(A(_, 2)) + 1, n = max(A(_, 3)) + 1;
// 		double ****B = ThreeD(m, l, k, n);

// 		for (int i = 0; i < A.nrow(); i++)
// 		{
// 			B[(int)A(i, 0)][(int)A(i, 1)][(int)A(i, 2)][A(i, 3)]++;
// 		}

// 		return unique(A);
// 		// return B[0][0][0];

// 	}  else if(A.ncol() == 5) {
// 		int l = max(A(_, 1)) + 1, k = max(A(_, 2)) + 1 n = max(A(_, 3)) + 1, o = max(A(_, 4)) + 1;
// 		double *****B = ThreeD(m, l, k, n, o);

// 		for (int i = 0; i < A.nrow(); i++)
// 		{
// 			B[(int)A(i, 0)][(int)A(i, 1)][(int)A(i, 2)][A(i, 3)][A(i, 4)]++;
// 		}

// 		return unique(A);
// 		// return B[0][0][0];

// 	} else {
// 		return unique(A);
// 	}
// }

// [[Rcpp::export]]
SEXP MySvalue(SEXP mat) {
	NumericMatrix A(mat);
	int hDim = A.ncol(), vDim = A.nrow();
	NumericVector sum(1, 0.0), pvalue(1, 0.0), out(2, 0.0);


	if (hDim == 2) {
		int l = max(A(_, 1)) + 1, m = max(A(_, 0)) + 1;
		double *x = OneD(m), *y = OneD(l), **z = TwoD(m, l);

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
				if (x[i] == 0 || y[j] == 0 || z[i][j] == 0)
					continue;

				sum[0] += 2.0 * z[i][j] * log( (z[i][j] * vDim) / (x[i] * y[j]) );
				// cout << z[i][j] << " * 2 * log( ( " << z[i][j] << " * " << vDim << " ) / ( " << x[i] << " * " << y[j] << " ) ) = " << 2 * z[i][j] * log( (z[i][j] * vDim) / (x[i] * y[j]) ) << endl;
			}
		}

		int DF = (m-1) * (l-1);
		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;

	} else if (hDim == 3) {
		int k = max(A(_, 2)) + 1, l = max(A(_, 1)) + 1, m = max(A(_, 0)) + 1;
		double *v = OneD(k), **x = TwoD(m, k), **y = TwoD(l, m), ***z = ThreeD(m, l, k);

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
					// cout << z[i][j][h] << " * 2 * log( ( " << z[i][j][h] << " * " << v[h] << " ) / ( " << x[i][h] << " * " << y[j][h] << " ) ) = " << 2 * z[i][j][h] * log( (z[i][j][h] * v[h]) / (x[i][h] * y[j][h]) ) << endl;
					if (x[i][h] == 0 || y[j][h] == 0 || z[i][j][h] == 0)
						continue;

					sum[0] += 2.0 * z[i][j][h] * log( (z[i][j][h] * v[h]) / (x[i][h] * y[j][h]) );
					// cout << sum[0] << endl;
				}
			}
		}

		// sum[0] *= 2;
		// int DF = Df(df);
		int DF = (m-1) * (l-1) * k;
		// cout << DF << endl;
		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;

	} else if (hDim == 4) {
		int n = max(A(_, 3)) + 1, k = max(A(_, 2)) + 1, l = max(A(_, 1)) + 1, m = max(A(_, 0)) + 1;
		double **v = TwoD(k, n), ***x = ThreeD(m, k, n), ***y = ThreeD(l, k, n), ****z = FourD(m, l, k, n);

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
						if (x[i][h][f] == 0 || y[j][h][f] == 0 || z[i][j][h][f] == 0)
							continue;
						sum[0] += 2.0 * z[i][j][h][f] * log( (z[i][j][h][f] * v[h][f]) / (x[i][h][f] * y[j][h][f]) );
						// cout << z[i][j][h][f] << " * 2 * log( ( " << z[i][j][h][f] << " * " << v[h][f] << " ) / ( " << x[i][h][f] << " * " << y[j][h][f] << " ) ) = " << 2 * z[i][j][h][f] * log( (z[i][j][h][f] * v[h][f]) / (x[i][h][f] * y[j][h][f]) ) << endl;
					}
				}
			}
		}

		int DF = (m-1) * (l-1) * k * (n-1);
		// cout << DF << endl;
		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;

	} else if (hDim == 5) {
		int o = max(A(_, 4)) + 1, n = max(A(_, 3)) + 1, k = max(A(_, 2)) + 1, l = max(A(_, 1)) + 1, m = max(A(_, 0)) + 1;
		double ***v = ThreeD(k, n, o), ****x = FourD(m, k, n, o), ****y = FourD(l, k, n, o), *****z = FiveD(m, l, k, n, o);

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
							if (x[i][h][f][e] == 0 || y[j][h][f][e] == 0 || z[i][j][h][f][e] == 0)
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
void Do(SEXP x) {
	NumericMatrix A(x);
	int k = max(A(_,0));
	NumericMatrix B(k, k);
	cout << k << endl;
}

// [[Rcpp::export]]
int Df(SEXP x) {
	double *DfSet = REAL(x);
	int df = (DfSet[0]-1)*(DfSet[1]-1);
	#pragma omp parallel
	{
		#pragma omp for
		{
			for (int i = 2; i < XLENGTH(x); i++)
			{
				df*=DfSet[i];
			}
		}
	}
	return df;
}

// [[Rcpp::export]]
NumericVector Statistics(SEXP x, SEXP y, SEXP z) {
	int n, m;
	// int *X = INTEGER(x);
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
		#pragma omp for collapse(2)
		{
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
				// cout << count[0] << " * 2 * log( ( " << count[0] << " * " << count[3] << " ) / ( " << count[1] << " * " << count[2] << " ) ) = " << 2 * count[0] * log( (count[0] * count[3]) / (count[1] * count[2]) ) << endl;
				sum[0] += 2 * count[0] * log( ( count[0] * count[3] ) / ( count[1] * count[2] ) );
				// cout << sum[0] << endl;
			}
		}
	}

	// df = Cardinality(A);
	int DF = Df(df);
	// cout << DF << endl;
	pvalue = pchisq(sum, DF, FALSE);
	out[0] = pvalue[0];
	out[1] = sum[0];

	return out;
}