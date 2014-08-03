#include <Rcpp.h>
#include <tr1/unordered_map>
#include <map>
#include <algorithm>
#include <time.h>

using namespace std;
using namespace std::tr1;
using namespace Rcpp;

int Hash(NumericVector& array, int i, bool skip) {
	int out = 0;
	for (i; i < array.size(); i++) {
		if (skip && i == 1)
			continue;
		out = 10*out + array[i];
	}
	
	return out;
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
NumericVector Svalue(NumericMatrix& A, const NumericVector& cardinality) {
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
		NumericMatrix B = T(A, hDim, vDim);
		unordered_map<int, int> Count;
		unordered_map<int, int> ReMap;
		int teta = 0;
		NumericVector tmp(hDim);

		for (int i = 0; i < vDim; i++)
		{
			tmp = B(_, i);
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

class Properties
{
public:
	Properties(SEXP);
	~Properties();

	template <typename T, int RTYPE> int colCardinality(const Vector<RTYPE>& x, unordered_map<T, int>& y);
	NumericVector Cardinality(SEXP);
	int Hash(NumericVector&, int);
	NumericVector It(NumericVector&, int);
	NumericMatrix partialMatrix(const NumericMatrix&, NumericVector&, int);
	NumericVector CorrespondingCardinality(NumericVector&, const NumericVector&, int);
	void Powerset(SEXP);
	
	int n;
	NumericVector cardinality;
	map<int, short> OK;
	unordered_map<int, double> P;
	unordered_map<int, double> S;
	List sums;
};

Properties::Properties(SEXP z) {
	cardinality = Cardinality(z);
	n = cardinality.size();
	Powerset(z);
}

Properties::~Properties() {}

template <typename T, int RTYPE> int Properties::colCardinality(const Vector<RTYPE>& x, unordered_map<T, int>& y) {
	int m = x.size();
	y.clear();
	for (int i = 0; i < m; i++)
		y[x[i]] = 1;

	return y.size();
}

NumericVector Properties::Cardinality(SEXP z) {
	NumericVector card;
	int i;
	switch (TYPEOF(z)) {
		case INTSXP:
		{
			IntegerMatrix A(z);
			IntegerVector x;
			unordered_map<int, int> y;
			for (i = 0; i < A.ncol(); i++) {
				x = A(_, i);
				card.push_back(colCardinality<int, INTSXP>(x, y));
			}
			break;
		}
		case REALSXP:
		{
			NumericMatrix A(z);
			NumericVector x;
			unordered_map<double, int> y;
			for (i = 0; i < A.ncol(); i++) {
				x = A(_, i);
				card.push_back(colCardinality<double, REALSXP>(x, y));
			}
			break;
		}
		default:
		{
			card.push_back(0);
			cout << "Error 02: Wrong datatype in matrix!" << endl;
		}
	}
	return card;
}

int Properties::Hash(NumericVector& array, int size) {
	int out = 0;
	for (int i = 0; i < size; i++)
		out = 10*out + (array[i] + 1);

	return out;
}

NumericVector Properties::It(NumericVector& array, int size) {
	NumericVector out;
	for (int i = 1; i <= size; i++)
		out.push_back(array[i]);

	return out;
}

NumericMatrix Properties::partialMatrix(const NumericMatrix& A, NumericVector& pa, int size) {
	NumericMatrix partMat(A.nrow(), size);

	for (int i = 0; i < size; i++)
		partMat(_, i) = A(_, pa(i));

	return partMat;
}

NumericVector Properties::CorrespondingCardinality(NumericVector& pa, const NumericVector& cardinality, int size) {
	NumericVector tmpCardinality;

	for (int i = 0; i < size; i++)
		tmpCardinality.push_back(cardinality[pa[i]]);

	return tmpCardinality;
}

// do {
// 		for (int i = 0; i < perm.size(); i++)
// 			cout << perm[i];
// 		cout << endl;
// 	} while (next_permutation(perm.begin(), perm.end()));

void Properties::Powerset (SEXP a) {
	NumericMatrix A(a), B;
	int n = A.ncol(), k, hash;
	NumericVector permutations(n), newCard, results, Iterator;
	for (int i = 0; i < n; i++)
		permutations[i] = i;

	do {
		for (k = 2; k <= n; k++) {
			Iterator = It(permutations, k);
			hash = Hash(Iterator, k);
			if (OK[hash] == 0) {
				B = partialMatrix(A, Iterator, k);
				newCard = CorrespondingCardinality(Iterator, cardinality, k);
				results = Svalue(B, newCard);

				if (results[0] < 0.05) {
					P[hash] = results[0];
					S[hash] = results[1];
					OK[hash] = 1;
				} else {
					OK[hash] = -1;
				}
			} else if (OK[hash] == 1) {
				continue;
			} else if (OK[hash] == -1) {
				break;
			}
		}
	} while (next_permutation(permutations.begin(), permutations.end()));


	// NumericVector stack(n), newCard, Iterator, output;
	// stack[0] = 0;

	// map<int, double> test;
	// while(1) {
	// 	if (stack[k] < n) {
	// 		stack[k + 1] = stack[k] + 1;
	// 		k++;
	// 	} else {
	// 		stack[k - 1]++;
	// 		k--;
	// 	}

	// 	if (k == 0)
	// 		break;

	// 	if (k > 1) {
 //  			do {
	// 			Iterator = It(stack, k);
	// 			hash = Hash(Iterator, k);
	// 			B = partialMatrix(A, Iterator, k);
	// 			newCard = CorrespondingCardinality(Iterator, cardinality, k);
	// 			output = Svalue(B, newCard);
	// 			P[hash] = output[0];
	// 			test[hash] = output[0];
	// 			S[hash] = output[1];
	// 		} while (next_permutation(Iterator.begin(), Iterator.end()));
	// 	}
	// }
	sums = List::create(wrap(P), wrap(S));
}

RCPP_MODULE(props) {
	class_<Properties>("Properties")
	.constructor<SEXP>()
	.field("dim", &Properties::n)
	.field("cardinality", &Properties::cardinality)
	.field("sums", &Properties::sums)
	.method("Cardinality", &Properties::Cardinality)
	.method("Powerset", &Properties::Powerset)
;}

int SetCols(NumericVector& cpc, int T, int X) {
	NumericVector pa;
	pa.push_back(T);
	pa.push_back(X);
	for (int i = 0; i < cpc.size(); i++)
		pa.push_back(cpc[i]);

	int hash = Hash(pa, 0, FALSE);

	return hash;
}

void UpdateCPC(List& CPC, int selected) {
	NumericVector tmp;
	if (CPC.size() == 0) {
		CPC.push_back(1);
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

void MaxMinHeuristic(Properties& A, NumericVector& cpc, List& CPCprops, NumericVector& variablesToTest, int T, double alpha) {//, double selectedBefore = -1.0, double minimum = 1.0) {
	NumericVector pvalue, rejectedInLastStep = as<NumericVector>(CPCprops[1]), temporaryMinimum = rejectedInLastStep, accepted;
	int pa;

	// for (int i = 0; i < variablesToTest.size(); i++)
	for (NumericVector::iterator it = variablesToTest.begin(); it != variablesToTest.end(); it++)
	{
		pa = SetCols(cpc, T, *it);
		if (A.OK[pa] == 1) {
			pvalue = NumericVector::create(A.P[pa], A.S[pa]);

			if (pvalue[0] < alpha) {
				if (pvalue[0] < rejectedInLastStep[1]) {
					temporaryMinimum = rejectedInLastStep;
					rejectedInLastStep = NumericVector::create(*it, pvalue[0], pvalue[1]);
				} else if (pvalue[0] == rejectedInLastStep[1] && pvalue[1] > rejectedInLastStep[2]) {
					rejectedInLastStep = NumericVector::create(*it, rejectedInLastStep[1], pvalue[1]);
				}
			} else {
				accepted.push_back(*it);
			}
		} else {
			continue;
		}
	}

	if (accepted.size() > 0)
		CPCprops[0] = accepted;
	if (rejectedInLastStep[0] != -1)
		CPCprops[1] = rejectedInLastStep;
	if (temporaryMinimum[0] != -1 && temporaryMinimum[0] != rejectedInLastStep[0])
		CPCprops[2] = temporaryMinimum;
}

// List Forward(NumericMatrix& A, int T, const NumericVector& cardinality, double alpha) {
// [[Rcpp::export]]
List Forward(Properties& A, int T, SEXP a) {
	NumericVector cpc, variablesToTest;

	int n = A.n;
	double *beta = REAL(a);
	double alpha = *beta;
	List CPC, CPCprops;

	NumericVector reject = NumericVector::create(-1.0, 1.0, 1.0), accepted = NumericVector::create(T), tmpAccepted;
	CPCprops.push_back(T);
	CPCprops.push_back(reject);
	CPCprops.push_back(reject);

	for (int i = 1; i <= n; i++)
		if (i != T)
			variablesToTest.push_back(i);

	UpdateCPC(CPC, 0);
	MaxMinHeuristic(A, cpc, CPCprops, variablesToTest, T, alpha);
	
	reject = as<NumericVector>(CPCprops[1]);
	tmpAccepted = as<NumericVector>(CPCprops[0]);

	if (reject[0] != -1.0) {
		UpdateCPC(CPC, reject[0]);
		cpc.push_back(reject[0]);
		accepted.push_back(reject[0]);
	}

	// for (int i = 0; i < tmpAccepted.size(); i++)
	for (NumericVector::iterator it = tmpAccepted.begin(); it != tmpAccepted.end(); it++)
		accepted.push_back(*it);

	variablesToTest = NumericVector::create();
	for (int i = 1; i <= n; i++)
	{
		if (!IsIn(accepted, i)) {
			variablesToTest.push_back(i);
		}
	}

	// if (variablesToTest.size() > 0) {
	while (variablesToTest.size() > 0) {
		CPCprops[1] = NumericVector::create(-1.0, 1.0, 1.0);
		MaxMinHeuristic(A, cpc, CPCprops, variablesToTest, T, alpha);

		reject = as<NumericVector>(CPCprops[1]);
		tmpAccepted = as<NumericVector>(CPCprops[0]);

		if (reject[0] != -1.0) {
			UpdateCPC(CPC, reject[0]);
			cpc.push_back(reject[0]);
			accepted.push_back(reject[0]);

			for (int i = 0; i < tmpAccepted.size(); i++)
				accepted.push_back(tmpAccepted[i]);

			variablesToTest = NumericVector::create();
			for (int i = 1; i <= n; i++)
			{
				if (!IsIn(accepted, i)) {
					variablesToTest.push_back(i);
				}
			}
		} else {
			variablesToTest = NumericVector::create();
		}
	}

	return CPC;
}

NumericVector Backward(Properties& A, SEXP a, List& CPC, int T) {
	double *alpha = REAL(a);
	NumericVector cpc = as<NumericVector>(CPC[CPC.size()-1]), rm, tmpCardinality, pvalue, fromCPC;
	int k = -1, pa;

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
				if (A.OK[pa] == 1) {
					pvalue = NumericVector::create(A.P[pa], A.S[pa]);

					if (pvalue[0] > *alpha && pvalue[0] != 1.0) {
						k = i;
						break;
					}
				} else {
					continue;
				}
			}
		}

		if (k != -1)
			cpc.erase(k);

		return cpc;
	}
}






// int getR(NumericVector& f) {
// 	unordered_map<int, int> Map;

// 	for (NumericVector::iterator it = f.begin(); it != f.end(); it++) {
// 		Map[*it] = 1;
// 	}

// 	return Map.size();
// }

// string int_array_to_string(NumericVector& array) { //  int int_array[], int size_of_array){
// 	ostringstream oss("");
// 	for (int temp = 0; temp < array.size(); temp++)
// 		oss << array[temp];
	
// 	return oss.str();
// }


// unordered_map<int, string> UniqueMap(SEXP x) {
// 	NumericMatrix A(x);
// 	NumericVector u;
// 	string key;
// 	unordered_map<string, int> Map;
// 	unordered_map<int, string> Uni;
// 	int n = Map.size(), k = 0;

// 	for (int i = 0; i < A.nrow(); i++)
// 	{
// 		u = A(i, _);
// 		key = int_array_to_string(u);
// 		Map[key] = 1;
// 		if (Map.size() != n) {
// 			Uni[k] = key;
// 			n++;
// 			k++;
// 		}
// 	}

// 	return Uni;
// }

// int getSingleN_ijk(double *vec, int dim, int k) {
// 	int count = 0;

// 	for (int i = 0; i < dim; i++)
// 	{
// 		if (vec[i] == k)
// 			count++;
// 	}

// 	return count;
// }

// int getVecN_ijk(double *vec, double *parentVec, int dim, int j, int k) {
// 	int count = 0;

// 	for (int i = 0; i < dim; i++)
// 	{
// 		if (vec[i] == k && parentVec[i] == j)
// 			count++;
// 	}

// 	return count;
// }

// int getMapN_ijk(double *vec, NumericMatrix& parentMatrix, unordered_map<int, string> parentMap, int dim, int j, int k) {
// 	int count = 0;
// 	NumericVector u;

// 	for (int i = 0; i < dim; i++)
// 	{
// 		u = parentMatrix(i, _);
// 		if (vec[i] == k && parentMap[j-1] == int_array_to_string(u))
// 			count++;
// 	}

// 	return count;
// }


// double ScoreNodeWithNoneParents(SEXP column, SEXP N, int r, double eta) {
// 	int *n = INTEGER(N), dim = XLENGTH(column);
// 	double *vec = REAL(column);
// 	double gammaJ = eta / 1.0, gammaK = eta / (1.0 * r);
// 	double rScore = 0.0, qScore = 0.0;
// 	int n_ijk = 0, n_ij = 0;

// 	for (int k = 1; k <= r; k++)
// 	{
// 		n_ijk = getSingleN_ijk(vec, dim, k);
// 		n_ij += n_ijk;
// 		if (n_ijk != 0)
// 			rScore += lgamma( n_ijk + gammaK ) - lgamma( gammaK );
// 	}

// 	qScore += lgamma( gammaJ ) - lgamma( n_ij + gammaJ ) + rScore;


// 	return qScore;
// }


// double ScoreNodeWithOneParent(SEXP Xi, SEXP Pa, SEXP N, int r, int q, double eta) {
// 	int *n = INTEGER(N), dim = XLENGTH(Xi);
// 	double *vec = REAL(Xi), *parentVec = REAL(Pa);
// 	double gammaJ = eta / q, gammaK = eta / (q * r);
// 	double rScore = 0.0, qScore = 0.0;
// 	int n_ij, n_ijk = 0;

// 	for (int j = 1; j <= q; j++)
// 	{
// 		n_ij = 0;
// 		rScore = 0.0;

// 		for (int k = 1; k <= r; k++)
// 		{
// 			n_ijk = getVecN_ijk(vec, parentVec, dim, j, k);
// 			n_ij += n_ijk;
// 			if (n_ijk != 0)
// 				rScore += lgamma( n_ijk + gammaK ) - lgamma( gammaK );
// 		}

// 		qScore += lgamma( gammaJ ) - lgamma( n_ij + gammaJ ) + rScore;
// 	}


// 	return qScore;
// }


// double ScoreNodeWithMoreParents(SEXP Xi, SEXP Pa, SEXP N, int r, int q, double eta) {
// 	int *n = INTEGER(N), dim = XLENGTH(Xi);
// 	double *vec = REAL(Xi);
// 	NumericMatrix Parent(Pa);
// 	unordered_map<int, string> parentMap = UniqueMap(Parent);
// 	double gammaJ = eta / q, gammaK = eta / ( q * r);
// 	double rScore = 0.0, qScore = 0.0;
// 	int n_ij, n_ijk;
// 	NumericVector N_ijk(r);

// 	for (int j = 1; j <= q; j++)
// 	{
// 		n_ij = 0;
// 		rScore = 0.0;

// 		for (int k = 1; k <= r; k++)
// 		{
// 			n_ijk = getMapN_ijk(vec, Parent, parentMap, dim, j, k);
// 			n_ij += n_ijk;
// 			if (n_ijk != 0)
// 				rScore += lgamma( n_ijk + gammaK ) - lgamma( gammaK );
// 		}

// 		qScore += lgamma( gammaJ ) - lgamma( n_ij + gammaJ ) + rScore;
// 	}


// 	return qScore;
// }

// NumericVector ReturnParents(int i, NumericMatrix& AdjMat, SEXP N) {
// 	NumericVector parents;
// 	int *n = INTEGER(N);
// 	for (int j = 0; j < *n; j++)
// 	{
// 		if (AdjMat(j, i) == 1) {
// 			parents.push_back(j);
// 		}
// 	}
// 	return parents;
// }

// NumericVector InitScore(NumericMatrix& A, SEXP N, NumericVector& R, double eta) {
// 	int *n = INTEGER(N);
// 	NumericVector g, score(*n);

// 	for (int i = 0; i < *n; i++) {
// 		g = A(_, i);
// 		score[i] = ScoreNodeWithNoneParents(g, N, R[i], eta);
// 	}

// 	return score;
// }


// NumericMatrix partialMatrix(const NumericMatrix& A, NumericVector& pa) {
// 	NumericMatrix partMat(A.nrow(), pa.size());

// 	for (int i = 0; i < pa.size(); i++)
// 	{
// 		partMat(_, i) = A(_, pa(i));
// 	}

// 	return partMat;
// }

// NumericVector ScoreGraph(const NumericMatrix& A, NumericMatrix& AdjMat, SEXP N, const NumericVector& R, double eta) {
// 	NumericVector childVector, parentVector, allParents, score;
// 	NumericMatrix parentMatrix;
// 	int q;

// 	for (int i = 0; i < AdjMat.ncol(); i++)
// 	{
// 		q = 1;
// 		allParents = ReturnParents(i, AdjMat, N);
// 		childVector = A(_, i);

// 		for (int k = 0; k < allParents.size(); k++)
// 			q *= R[allParents[k]];

// 		if (allParents.size() == 0) {
// 			score.push_back(ScoreNodeWithNoneParents(childVector, N, R[i], eta));
// 		} else if (allParents.size() == 1) {
// 			parentVector = A(_, allParents[0]);
// 			score.push_back(ScoreNodeWithOneParent(childVector, parentVector, N, R[i], q, eta));
// 		} else {
// 			parentMatrix = partialMatrix(A, allParents);
// 			score.push_back(ScoreNodeWithMoreParents(childVector, parentMatrix, N, R[i], q, eta));
// 		}
// 	}

// 	return score;
// }

// void SettingEdges(const NumericMatrix& A, NumericMatrix& AdjMat, NumericVector& scores, const List& PC, SEXP N, const NumericVector& R, double eta) {
// 	NumericVector childVector, parentVector, allParents, tmp;
// 	NumericMatrix parentMatrix;
// 	double add, before, after;
// 	int q;

// 	for (int i = 0; i < PC.size(); i++)
// 	{
// 		NumericVector pc = as<NumericVector>(PC[i]);
// 		// q = 1;
// 		// add = scores[i];

// 		for (int j = 0; j < pc.size(); j++) {

// 			if (AdjMat(i, pc[j] - 1) == 0 && AdjMat(pc[j] - 1, i) == 0) {
				
// 				AdjMat(i, pc[j] - 1) = 1;

// 				before = sum(scores);
// 				tmp = ScoreGraph(A, AdjMat, N, R, eta);
// 				after = sum(tmp);

// 				if (after > before) {
// 					scores = tmp;
// 				} else {
// 					AdjMat(i, pc[j] - 1) = 0;
// 				}
// 			}
// 		}
// 	}
// }

// void AddReverseDelete(const NumericMatrix& A, NumericMatrix& AdjMat, NumericVector& scores, const List& PC, SEXP N, const NumericVector& R, double eta) {
// 	srand(time(NULL));
// 	NumericVector tmp;
// 	double before, after;
// 	int rnd;

// 	for (int i = 0; i < PC.size(); i++)
// 	{
// 		NumericVector pc = as<NumericVector>(PC[i]);

// 		for (int j = 0; j < pc.size(); j++) {

// 			rnd = (rand() % 10) + 1;

// 			if (AdjMat(i, pc[j] - 1) == 0 && AdjMat(pc[j] - 1, i) == 0) {
				
// 				AdjMat(pc[j] - 1, i) = 1;

// 				before = sum(scores);
// 				tmp = ScoreGraph(A, AdjMat, N, R, eta);
// 				after = sum(tmp);

// 				if (after > before) {
// 					scores = tmp;
// 				} else {
// 					AdjMat(pc[j] - 1, i) = 0;
// 				}
// 			} else if (AdjMat(i, pc[j] - 1) == 1 && rnd > 7) {
				
// 				AdjMat(i, pc[j] - 1) = 0;

// 				before = sum(scores);
// 				tmp = ScoreGraph(A, AdjMat, N, R, eta);
// 				after = sum(tmp);

// 				if (after > before) {
// 					scores = tmp;
// 				} else {
// 					AdjMat(i, pc[j] - 1) = 1;
// 				}
// 			}  else if (AdjMat(pc[j] - 1, i) == 1 && rnd <= 7) {

// 				AdjMat(pc[j] - 1, i) = 0;
// 				AdjMat(i, pc[j] - 1) = 1;

// 				before = sum(scores);
// 				tmp = ScoreGraph(A, AdjMat, N, R, eta);
// 				after = sum(tmp);

// 				if (after > before) {
// 					scores = tmp;
// 				} else {
// 					AdjMat(pc[j] - 1, i) = 1;
// 					AdjMat(i, pc[j] - 1) = 0;
// 				}
// 			}
// 		}
// 	}
// }

// SEXP BDeu(SEXP x, SEXP y, SEXP z, SEXP N) {
// 	int *n= INTEGER(N);
// 	double eta = 2.2;
// 	NumericVector R(y), scores, tmpScores;
// 	NumericMatrix A(x), AdjMat(*n, *n), tmpAdjMat(*n, *n);
// 	List PC(z);

// 	scores = InitScore(A, N, R, eta);
// 	SettingEdges(A, AdjMat, scores, PC, N, R, eta);
// 	AddReverseDelete(A, tmpAdjMat, tmpScores, PC, N, R, eta);
	
// 	for (int i = 0; i < 5; i++)
// 	{
// 		tmpScores = scores;
// 		tmpAdjMat = AdjMat;
// 		AddReverseDelete(A, tmpAdjMat, tmpScores, PC, N, R, eta);
		
// 		if (sum(scores) >= sum(tmpScores)) {
// 			break;
// 		} else {
// 			AdjMat = tmpAdjMat;
// 			scores = tmpScores;
// 		}
// 	}

// 	cout << sum(scores) << endl;
	
// 	return AdjMat;
// }

// [[Rcpp::export]]
List MMHC(SEXP x, SEXP a) {
	NumericVector cpc, pc, tmppc;
	List CPC, PC, tmpPC;
	NumericVector count;
	Properties A(x);
	
	for (int T = 1; T <= 5; T++)
	{
		CPC = Forward(A, T, a);
		cpc = as<NumericVector>(CPC[CPC.size()-1]);
		if (cpc.size() == 0) {
			PC.push_back(R_NilValue);
		} else {
			pc = Backward(A, a, CPC, T);

			for (int j = 0; j < pc.size(); j++)
			{
				tmpPC = Forward(A, (int)pc[j], a);
				tmppc = Backward(A, a, tmpPC, (int)pc[j]);
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
			PC.push_back(pc);
		}
	}

	return PC;
}