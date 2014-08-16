#include <Rcpp.h>
#include <tr1/unordered_map>

using namespace std;
using namespace std::tr1;
using namespace Rcpp;

int Hash(NumericVector& array, int i, bool skip) {
	int out = 0;
	for (i; i < array.size(); i++)
		if (skip && i == 1)
			continue;
		out += 10*out + array[i];

	return out;
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
SEXP Svalue(NumericMatrix& A, const NumericVector& cardinality) {
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

	} else if (hDim == 4) {
		int n = cardinality[3], k = cardinality[2], l = cardinality[1], m = cardinality[0];
		double **v = TwoD(k, n), ***x = ThreeD(m, k, n), ***y = ThreeD(l, k, n), ****z = FourD(m, l, k, n);

		for (int i = 0; i < vDim; i++)
		{
			v[(int)A(i, 2) - 1][(int)A(i, 3) - 1]++;
			x[(int)A(i, 0) - 1][(int)A(i, 2) - 1][(int)A(i, 3) - 1]++;
			y[(int)A(i, 1) - 1][(int)A(i, 2) - 1][(int)A(i, 3) - 1]++;
			z[(int)A(i, 0) - 1][(int)A(i, 1) - 1][(int)A(i, 2) - 1][(int)A(i, 3) - 1]++;
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
			v[(int)A(i, 2) - 1][(int)A(i, 3) - 1][(int)A(i, 4) - 1]++;
			x[(int)A(i, 0) - 1][(int)A(i, 2) - 1][(int)A(i, 3) - 1][(int)A(i, 4) - 1]++;
			y[(int)A(i, 1) - 1][(int)A(i, 2) - 1][(int)A(i, 3) - 1][(int)A(i, 4) - 1]++;
			z[(int)A(i, 0) - 1][(int)A(i, 1) - 1][(int)A(i, 2) - 1][(int)A(i, 3) - 1][(int)A(i, 4) - 1]++;
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
		int acKey, bcKey, abcKey, cKey;
		unordered_map<int, int> Count;
		unordered_map<int, int> ReMap;
		int teta = 0;
		NumericVector tmp(hDim);

		for (int i = 0; i < hDim; i++)
		{
			tmp = A(i, _);
			abcKey = Hash(tmp, 0, FALSE);
			bcKey = Hash(tmp, 1, FALSE);
			cKey = Hash(tmp, 2, FALSE);
			acKey = Hash(tmp, 0, TRUE);
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

template <typename T, int RTYPE> int colCardinality(const Vector<RTYPE>& x, unordered_map<T, int>& y) {
	int m = x.size();
	y.clear();
	for (int i = 0; i < m; i++)
		y[x[i]] = 1;

	return y.size();
}

void Cardinality(NumericMatrix& A, NumericVector& cardinality) {
	NumericVector x;
	unordered_map<double, int> y;
	for (int i = 0; i < A.ncol(); i++)
	{
		x = A(_, i);
		cardinality.push_back(colCardinality<double, REALSXP>(x, y));
	}
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
		if (k == cpc.size() && i != *T)
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

NumericVector SetCols(NumericVector& cpc, int T, int X) {
	NumericVector pa;
	pa.push_back(X);
	pa.push_back(T);
	for (int i = 0; i < cpc.size(); i++)
	{
		pa.push_back(cpc[i]);
	}
	return pa;
}

void UpdateCPC(List& CPC, int selected) {
	NumericVector tmp;
	if (CPC.size() == 0) {
		CPC.push_back(0);
		CPC.push_back(R_NilValue);
	} else if (CPC.size() == 2) {
		CPC[0] = 2;
		CPC.push_back(selected);
	} else {
		int CPCLength = CPC.size();
		CPC[0] = CPCLength;
		CPC.push_back(selected);
		for (int i = 2; i < CPCLength; i++)
		{
			tmp = as<NumericVector>(CPC[i]);
			tmp.push_back(selected);
			CPC.push_back(tmp);
		}
	}
}

bool IsIn(NumericVector& x, double y) {
	bool out = FALSE;
	for (int i = 0; i < x.size(); i++)
	{
		if (x[i] == y) {
			out = TRUE;
			break;
		}
	}

	return out;
}

void MaxMinHeuristic(NumericMatrix& A, NumericVector& cpc, List& CPCprops, NumericVector& variablesToTest, const NumericVector& cardinality, int T, double alpha) {//, double selectedBefore = -1.0, double minimum = 1.0) {
	NumericMatrix statisticMatrix;
	NumericVector pa, pvalue, tmpCardinality, rejectedInLastStep = as<NumericVector>(CPCprops[1]), temporaryMinimum = rejectedInLastStep, accepted;

	for (int i = 0; i < variablesToTest.size(); i++)
	{
		pa = SetCols(cpc, T, variablesToTest[i]);
		tmpCardinality = CorrespondingCardinality(pa, cardinality);
		statisticMatrix = partialMatrix(A, pa);
		pvalue = Svalue(statisticMatrix, tmpCardinality);

		if (pvalue[0] < alpha) {
			if (pvalue[0] < rejectedInLastStep[1]) {
				temporaryMinimum = rejectedInLastStep;
				rejectedInLastStep = NumericVector::create(variablesToTest[i], pvalue[0], pvalue[1]);
			} else if (pvalue[0] == rejectedInLastStep[1] && pvalue[1] > rejectedInLastStep[2]) {
				rejectedInLastStep = NumericVector::create(variablesToTest[i], rejectedInLastStep[1], pvalue[1]);
			}
		} else {
			accepted.push_back(variablesToTest[i]);
		}
	}

	if (accepted.size() > 0)
		CPCprops[0] = accepted;
	if (rejectedInLastStep[0] != -1)
		CPCprops[1] = rejectedInLastStep;
	if (temporaryMinimum[0] != -1 && temporaryMinimum[0] != rejectedInLastStep[0])
		CPCprops[2] = temporaryMinimum;
}

void CompatibilityToR(NumericVector& cpc) {
	for (int i = 0; i < cpc.size(); i++)
		cpc[i]++;
}

// List Forward(NumericMatrix& A, int T, const NumericVector& cardinality, double alpha) {
List Forward(NumericMatrix& A, int T, NumericVector& cardinality, SEXP a) {
	// NumericMatrix A(x);
	NumericVector cpc, variablesToTest;//, cardinality(z);
	// int *t = INTEGER(y), n = A.ncol();
	// int T = *t - 1;
	int n = A.ncol();
	double *beta = REAL(a);
	double alpha = *beta;
	List CPC, CPCprops;
	// *T = *T - 1;

	NumericVector reject = NumericVector::create(-1.0, 1.0, 1.0), accepted = NumericVector::create(T), tmpAccepted;
	CPCprops.push_back(T);
	CPCprops.push_back(reject);
	CPCprops.push_back(reject);

	for (int i = 0; i < n; i++)
		if (i != T)
			variablesToTest.push_back(i);

	UpdateCPC(CPC, 0);
	MaxMinHeuristic(A, cpc, CPCprops, variablesToTest, cardinality, T, alpha);
	
	reject = as<NumericVector>(CPCprops[1]);
	tmpAccepted = as<NumericVector>(CPCprops[0]);

	if (reject[0] != -1.0) {
		UpdateCPC(CPC, reject[0]);
		cpc.push_back(reject[0]);
		accepted.push_back(reject[0]);
	}

	for (int i = 0; i < tmpAccepted.size(); i++)
		accepted.push_back(tmpAccepted[i]);

	variablesToTest = NumericVector::create();
	for (int i = 0; i < n; i++)
	{
		if (!IsIn(accepted, i)) {
			variablesToTest.push_back(i);
		}
	}

	// if (variablesToTest.size() > 0) {
	while (variablesToTest.size() > 0) {
		CPCprops[1] = NumericVector::create(-1.0, 1.0, 1.0);
		MaxMinHeuristic(A, cpc, CPCprops, variablesToTest, cardinality, T, alpha);

		reject = as<NumericVector>(CPCprops[1]);
		tmpAccepted = as<NumericVector>(CPCprops[0]);

		if (reject[0] != -1.0) {
			UpdateCPC(CPC, reject[0]);
			cpc.push_back(reject[0]);
			accepted.push_back(reject[0]);

			for (int i = 0; i < tmpAccepted.size(); i++)
				accepted.push_back(tmpAccepted[i]);

			variablesToTest = NumericVector::create();
			for (int i = 0; i < n; i++)
			{
				if (!IsIn(accepted, i)) {
					variablesToTest.push_back(i);
				}
			}
		} else {
			variablesToTest = NumericVector::create();
		}
	}

	// CompatibilityToR(cpc);

	return CPC;
}

NumericVector Backward(NumericMatrix& A, NumericVector& cardinality, SEXP a, List& CPC, int T) {
	// NumericMatrix A(x);
	// NumericVector cardinality(z);
	double *alpha = REAL(a);
	NumericVector cpc = as<NumericVector>(CPC[CPC.size()-1]), rm, pa, tmpCardinality, pvalue, fromCPC;
	NumericMatrix statisticMatrix;
	int k = -1;

	if (cpc.size() == 1) {
		return cpc;
	} else {

		for (int i = 0; i < cpc.size(); i++)
		{
			for (int j = 1; j < CPC.size(); j++)
			{
				if (j == 1) {
					fromCPC = NumericVector::create();
				} else {
					fromCPC = as<NumericVector>(CPC[j]);
				}

				pa = SetCols(fromCPC, T, cpc[i]);
				tmpCardinality = CorrespondingCardinality(pa, cardinality);
				statisticMatrix = partialMatrix(A, pa);
				pvalue = Svalue(statisticMatrix, tmpCardinality);

				if (pvalue[0] > *alpha && pvalue[0] != 1.0) {
					k = i;
					break;
				}
			}
		}

		if (k != -1)
			cpc.erase(k);

		return cpc;
	}
}

// [[Rcpp::export]]
List MMPC(SEXP x, SEXP a) {
	NumericMatrix A(x);
	NumericVector cpc, pc, tmppc;
	List CPC, PC, tmpPC;
	NumericVector count, cardinality;
	Cardinality(A, cardinality);
	
	for (int T = 0; T < 5; T++)
	{
		CPC = Forward(A, T, cardinality, a);
		cpc = as<NumericVector>(CPC[CPC.size()-1]);
		if (cpc.size() == 0) {
			PC.push_back(R_NilValue);
		} else {
			pc = Backward(A, cardinality, a, CPC, T);

			for (int j = 0; j < pc.size(); j++)
			{
				tmpPC = Forward(A, (int)pc[j], cardinality, a);
				tmppc = Backward(A, cardinality, a, tmpPC, (int)pc[j]);
				if (!(IsIn(tmppc, T))) {
					int k = pc.size();
					for (int l = 0; l < k; l++)
					{
						if (pc[l] == pc[j]) {
							pc.erase(l);
						}
					}
				}
			}

			CompatibilityToR(pc);
			PC.push_back(pc);
		}

		// cpc = Backward(A, CPC, cardinality, T, 0.05);
		// PC.push_back(cpc);
	}

	// PC = Forward(x, y, z, a);
	// PC.push_back(cpc);

	return PC;
}