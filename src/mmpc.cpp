#include <omp.h>
#include <R.h>
#include <Rmath.h>
#include <vector>
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

string HashC(NumericVector& array, int i, bool skip) {
	ostringstream oss("");
	for (i; i < array.size(); i++) {
		if (skip && i == 1)
			continue;
		oss << array[i];
	}
	
	return oss.str();
}

NumericMatrix T(NumericMatrix& A, int row, int col) {
	NumericMatrix B(row, col);
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			B(i, j) = A(j, i);
	return B;
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

SEXP MySvalue(NumericMatrix& A, const NumericVector& cardinality) {
	int hDim = A.ncol(), vDim = A.nrow();
	NumericVector sum(1, 0.0), pvalue(1, 0.0), out(2, 0.0);

	if (hDim == 2) {
		int l = cardinality[1], m = cardinality[0];
		double *x = OneD(m), *y = OneD(l), **z = TwoD(m, l);

		for (int i = 0; i < vDim; i++)
		{
			x[(int)A(i, 0) - 1]++;
			y[(int)A(i, 1) - 1]++;
			z[(int)A(i, 0) - 1][(int)A(i, 1) - 1]++;
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
			v[(int)A(i, 2) - 1]++;
			x[(int)A(i, 0) - 1][(int)A(i, 2) - 1]++;
			y[(int)A(i, 1) - 1][(int)A(i, 2) - 1]++;
			z[(int)A(i, 0) - 1][(int)A(i, 1) - 1][(int)A(i, 2) - 1]++;
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

	} else {
		string acKey, bcKey, abcKey, cKey;
		NumericMatrix B = T(A, hDim, vDim);
		unordered_map<string, int> Count;
		unordered_map<int, string> ReMap;
		int teta = 0;
		NumericVector tmp(hDim);

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

void Cardinality(NumericMatrix& A, NumericVector& cardinality) {
	for (int i = 0; i < A.ncol(); i++)
		cardinality[i] = max(A(_, i));
}

NumericMatrix partialMatrix(const NumericMatrix& A, NumericVector& pa) {
	NumericMatrix partMat(A.nrow(), pa.size());

	for (int i = 0; i < pa.size(); i++)
		partMat(_, i) = A(_, pa(i));

	return partMat;
}

void ElementsToTest(NumericVector& variablesToTest, NumericVector& cpc, int *T, int n) {
	variablesToTest = NumericVector::create();
	for (int i = 0; i < n; i++)
	{
		int k = 0;
		for (int j = 0; j < cpc.size(); j++)
		{
			if (i != cpc[j])
				k++;
		}
		if (k == cpc.size() && i+1 != *T)
			variablesToTest.push_back(i);
	}
}

NumericVector CorrespondingCardinality(NumericVector& pa, const NumericVector& cardinality) {
	NumericVector tmpCardinality;

	for (int i = 0; i < pa.size(); i++)
	{
		tmpCardinality.push_back(cardinality[pa[i]]);
	}

	return tmpCardinality;
}

NumericVector SetCols(NumericVector& cpc, int *T, int X) {
	NumericVector pa;
	pa.push_back(*T-1);
	pa.push_back(X);
	for (int i = 0; i < cpc.size(); i++)
	{
		pa.push_back(cpc[i]);
	}
	return pa;
}

void MaxMinHeuristic(NumericMatrix& A, NumericVector& cpc, NumericVector& variablesToTest, NumericVector& temporaryMinimum, const NumericVector& cardinality, int *T, double *alpha, int n, double selectedBefore = -1.0, double minimum = 1.0) {
	NumericMatrix statisticMatrix;
	NumericVector pa, acceptedInThisStep, pvalue, tmpCardinality, variablesNotToTest, rejectedInLastStep = NumericVector::create(selectedBefore, minimum, 1.0);
	temporaryMinimum = rejectedInLastStep;

	for (int i = 0; i < variablesToTest.size(); i++)
	{
		pa = SetCols(cpc, T, variablesToTest[i]);
		tmpCardinality = CorrespondingCardinality(pa, cardinality);

		statisticMatrix = partialMatrix(A, pa);
		pvalue = MySvalue(statisticMatrix, tmpCardinality);

		if (pvalue[0] < *alpha) {
			if (pvalue[0] < rejectedInLastStep[1]) {
				temporaryMinimum = rejectedInLastStep;
				rejectedInLastStep = NumericVector::create(variablesToTest[i], pvalue[0], pvalue[1]);
			} else if (pvalue[0] == rejectedInLastStep[1] && pvalue[1] > rejectedInLastStep[2]) {
				rejectedInLastStep = NumericVector::create(variablesToTest[i], rejectedInLastStep[1], pvalue[1]);
			}
		} else {
			variablesNotToTest.push_back(variablesToTest[i]);
		}
	}

	if (rejectedInLastStep[0] != -1)
		cpc.push_back(rejectedInLastStep[0]);

	for (int i = 0; i < cpc.size(); i++)
	{
		cout << cpc[i] << endl;
	}
	cout << endl;

	for (int i = 0; i < cpc.size(); i++)
		variablesNotToTest.push_back(cpc[i]);

	ElementsToTest(variablesToTest, variablesNotToTest, T, n);

}

void CompatibilityToR(NumericVector& cpc) {
	for (int i = 0; i < cpc.size(); i++)
	{
		cpc[i]++;
	}
}

// [[Rcpp::export]]
NumericVector Forward(SEXP x, SEXP y, SEXP z, SEXP a) {
	NumericMatrix A(x), statisticMatrix;
	NumericVector cpc, variablesToTest, temporaryMinimum, cardinality(z);
	int *T = INTEGER(y), n = A.ncol();
	double *alpha = REAL(a);

	for (int i = 0; i < n; i++)
		if (i+1 != *T)
			variablesToTest.push_back(i);

	MaxMinHeuristic(A, cpc, variablesToTest, temporaryMinimum, cardinality, T, alpha, n);

	while (variablesToTest.size() != 0) {
		if (temporaryMinimum[0] == -1) {
			MaxMinHeuristic(A, cpc, variablesToTest, temporaryMinimum, cardinality, T, alpha, n);
		} else {
			MaxMinHeuristic(A, cpc, variablesToTest, temporaryMinimum, cardinality, T, alpha, n, temporaryMinimum[0], temporaryMinimum[1]);
		}
	}

	CompatibilityToR(cpc);

	return cpc;
}

void UpdateCPC(List& cpc, int selected) {
	NumericVector tmp;
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
SEXP initEmptyList(SEXP n) {
	int *m = INTEGER(n);
	List PC;

	for (int i = 0; i < *m; i++)
		PC.push_back(R_NilValue);

	return PC;
}

// [[Rcpp::export]]
int getR(SEXP x) {
	IntegerVector f(x);
	unordered_map<int, int> Map;

	for (IntegerVector::iterator it = f.begin(); it != f.end(); it++) {
		Map[*it] = 1;
	}

	return Map.size();
}

// [[Rcpp::export]]
bool Identical(SEXP x, SEXP y) {
	int *f = INTEGER(x), *g = INTEGER(y);
	int n = XLENGTH(x), m = XLENGTH(y), i = 0;
	bool id = TRUE;

	if (n != m) {
		id = FALSE;
	} else {
		while (f[i] == g[i]) {
			i++;
		}
	}
	if (i != n)
		id = FALSE;

	return id;
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
SEXP ToCrossOutC(SEXP x, SEXP y) {
	IntegerVector A(x), B(y), C;

	for (IntegerVector::iterator it = A.begin() ; it != A.end(); it++) {
		if (!(find(B.begin(), B.end(), *it) != B.end())) {
			C.push_back(*it);
		}
	}
	return C;
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