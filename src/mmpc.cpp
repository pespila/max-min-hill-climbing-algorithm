#include "mmhc.h"

MMHC::MMHC() {}

MMHC::MMHC(SEXP x) {
    DataFrame B(x);
	this->vDim = B.nrows();
	this->hDim = B.length();
	this->A = IntegerMatrix(this->vDim, this->hDim);
	this->graph = IntegerMatrix(this->hDim, this->hDim);
	this->scores = NumericVector(this->hDim);
	IntegerVector tmp;
	for (int i = 0; i < this->hDim; i++) {
		tmp = B[i];
		this->A(_, i) = tmp;
	}
	this->alpha = 0.05;
	this->eta = 1.0;
	// this->eta = (double)(sum(this->cardinality))/(double)(this->cardinality.size());
	this->cardinality = IntegerVector(this->hDim, 0);
	Cardinality();
    this->maxi = max(this->cardinality);
    this->D1 = this->OneD(this->maxi);
    this->D11 = this->OneD(this->maxi);
    this->D2 = this->TwoD(this->maxi, this->maxi);
    this->D22 = this->TwoD(this->maxi, this->maxi);
    this->D3 = this->ThreeD(this->maxi, this->maxi, this->maxi);
    this->D33 = this->ThreeD(this->maxi, this->maxi, this->maxi);
    this->D4 = this->FourD(this->maxi, this->maxi, this->maxi, this->maxi);
    this->D44 = this->FourD(this->maxi, this->maxi, this->maxi, this->maxi);
    this->D5 = this->FiveD(this->maxi, this->maxi, this->maxi, this->maxi, this->maxi);
}

MMHC::~MMHC() {
    this->maxi = max(this->cardinality);
    for (int i = 0; i < this->maxi; i++)
    {
    	for (int j = 0; j < this->maxi; j++)
		{
			for (int a = 0; a < this->maxi; a++)
			{
				for (int b = 0; b < this->maxi; b++)
				{
					free(D5[i][j][a][b]);
				}
				free(D5[i][j][a]);
			}
			free(D5[i][j]);
		}
		free(D5[i]);
	}
	free(D5);
    for (int i = 0; i < this->maxi; i++)
    {
        for (int j = 0; j < this->maxi; j++)
		{
			for (int a = 0; a < this->maxi; a++)
			{
				free(D4[i][j][a]);
			}
			free(D4[i][j]);
		}
		free(D4[i]);
	}
	free(D4);
    for (int i = 0; i < this->maxi; i++)
    {
    	for (int j = 0; j < this->maxi; j++)
		{
			for (int a = 0; a < this->maxi; a++)
			{
				free(D44[i][j][a]);
			}
			free(D44[i][j]);
		}
		free(D44[i]);
	}
	free(D44);
    for (int i = 0; i < this->maxi; i++)
    {
        for (int j = 0; j < this->maxi; j++)
		{
			free(D3[i][j]);
		}
		free(D3[i]);
	}
	free(D3);
    for (int i = 0; i < this->maxi; i++)
    {
    	for (int j = 0; j < this->maxi; j++)
		{
			free(D33[i][j]);
		}
		free(D33[i]);
	}
	free(D33);
    for (int i = 0; i < this->maxi; i++)
    	free(D2[i]);
	free(D2);
    for (int i = 0; i < this->maxi; i++)
		free(D22[i]);
	free(D22);
	free(D1);
	free(D11);
}

template <typename T, int RTYPE> int MMHC::colCardinality(const Vector<RTYPE>& x, unordered_map<T, int>& y) {
	int m = x.size();
	y.clear();
	for (int i = 0; i < m; i++)
		y[x[i]] = 1;

	return y.size();
}

void MMHC::Cardinality() {
	IntegerVector x;
	unordered_map<int, int> y;
	for (int i = 0; i < this->hDim; i++) {
		x = this->A(_, i);
		this->cardinality[i] = colCardinality<int, INTSXP>(x, y);
	}
}

int MMHC::Hash(IntegerVector& array, int i, bool skip) {
	int out = 0;
	for (i; i < array.size(); i++) {
		if (skip && i == 1)
			continue;
		out += 10*out + array[i];
	}

	return out;
}

double *MMHC::OneD(double x) {
	double *matrix = (double*)calloc(x, sizeof(double));
	for (int i = 0; i < x; i++)
		matrix[i] = 0;
	return matrix;
}

double **MMHC::TwoD(double x, double y) {
	double **matrix = (double**)calloc(x, sizeof(double*));
	for (int i = 0; i < x; i++)	{
		matrix[i] = (double*)calloc(y, sizeof(double));
		for (int j = 0; j < y; j++)
			matrix[i][j] = 0;
	}
	return matrix;
}

double ***MMHC::ThreeD(double x, double y, double z) {
	double ***matrix = (double***)calloc(x, sizeof(double*));
	for (int i = 0; i < x; i++)	{
		matrix[i] = (double**)calloc(y, sizeof(double*));
		for (int j = 0; j < y; j++)	{
			matrix[i][j] = (double*)calloc(z, sizeof(double));
			for (int k = 0; k < z; k++)
				matrix[i][j][k] = 0;
		}
	}
	return matrix;
}

double ****MMHC::FourD(double x, double y, double z, double a) {
	double ****matrix = (double****)calloc(x, sizeof(double*));
	for (int i = 0; i < x; i++)	{
		matrix[i] = (double***)calloc(y, sizeof(double*));
		for (int j = 0; j < y; j++)	{
			matrix[i][j] = (double**)calloc(z, sizeof(double*));
			for (int k = 0; k < z; k++)	{
				matrix[i][j][k] = (double*)calloc(a, sizeof(double));
				for (int l = 0; l < a; l++)
					matrix[i][j][k][l] = 0;
			}
		}
	}
	return matrix;
}

double *****MMHC::FiveD(double x, double y, double z, double a, double b) {
	double *****matrix = (double*****)calloc(x, sizeof(double*));
	for (int i = 0; i < x; i++)	{
		matrix[i] = (double****)calloc(y, sizeof(double*));
		for (int j = 0; j < y; j++)	{
			matrix[i][j] = (double***)calloc(z, sizeof(double*));
			for (int k = 0; k < z; k++)	{
				matrix[i][j][k] = (double**)calloc(a, sizeof(double*));
				for (int l = 0; l < a; l++)	{
					matrix[i][j][k][l] = (double*)calloc(b, sizeof(double));
					for (int m = 0; m < b; m++)
						matrix[i][j][k][l][m] = 0;
				}
			}
		}
	}
	return matrix;
}

NumericVector MMHC::Svalue(IntegerMatrix& A, const IntegerVector& cardinality) {
	int hDim = A.ncol(), vDim = A.nrow();
	NumericVector sum(1, 0.0), pvalue(1, 0.0), out(2, 0.0);

	if (hDim == 2) {
		int l = cardinality[1], m = cardinality[0];
//		double *x = OneD(m), *y = OneD(l), **z = TwoD(m, l);

		for (int i = 0; i < vDim; i++) {
			this->D1[A(i, 0) - 1]++;
			this->D11[A(i, 1) - 1]++;
			this->D2[A(i, 0) - 1][A(i, 1) - 1]++;
		}

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < l; j++) {
				if (this->D1[i] < 1 || this->D11[j] < 1 || this->D2[i][j] < 1)
					continue;
				sum[0] += 2.0 * this->D2[i][j] * log( (this->D2[i][j] * vDim) / (this->D1[i] * this->D11[j]) );
                this->D2[i][j] = 0;
			}
            this->D1[i] = 0;
		}
        
        for (int i = 0; i < l; i++)
            this->D11[i] = 0;

		int DF = cardinality[0] * cardinality[1];
		for (int i = 2; i < cardinality.size(); i++)
			DF *= (cardinality[i] + 1);

		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;

	}  else if (hDim == 3) {
		int k = cardinality[2], l = cardinality[1], m = cardinality[0];
//		double *v = OneD(k), **x = TwoD(m, k), **y = TwoD(l, k), ***z = ThreeD(m, l, k);

		for (int i = 0; i < vDim; i++) {
			this->D1[A(i, 2) - 1]++;
			this->D2[A(i, 0) - 1][A(i, 2) - 1]++;
			this->D22[A(i, 1) - 1][A(i, 2) - 1]++;
			this->D3[A(i, 0) - 1][A(i, 1) - 1][A(i, 2) - 1]++;
		}

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < l; j++) {
				for (int h = 0; h < k; h++) {
					if (this->D2[i][h] < 1 || this->D22[j][h] < 1 || this->D3[i][j][h] < 1 || this->D1[h] < 1)
						continue;
					sum[0] += 2.0 * this->D3[i][j][h] * log( (this->D3[i][j][h] * this->D1[h]) / (this->D2[i][h] * this->D22[j][h]) );
                    this->D3[i][j][h] = 0;
				}
			}
		}
		for (int i = 0; i < this->maxi; i++)
		{
			for (int j = 0; j < this->maxi; j++)
			{
				this->D2[i][j] = 0;
                this->D22[i][j] = 0;
			}
            this->D1[i] = 0;
		}

		int DF = cardinality[0] * cardinality[1];
		for (int i = 2; i < cardinality.size(); i++)
			DF *= (cardinality[i] + 1);

		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;

	} else if (hDim == 4) {
		int n = cardinality[3], k = cardinality[2], l = cardinality[1], m = cardinality[0];
//		double **v = TwoD(k, n), ***x = ThreeD(m, k, n), ***y = ThreeD(l, k, n), ****z = FourD(m, l, k, n);

		for (int i = 0; i < vDim; i++) {
			this->D2[A(i, 2) - 1][A(i, 3) - 1]++;
			this->D3[A(i, 0) - 1][A(i, 2) - 1][A(i, 3) - 1]++;
			this->D33[A(i, 1) - 1][A(i, 2) - 1][A(i, 3) - 1]++;
			this->D4[A(i, 0) - 1][A(i, 1) - 1][A(i, 2) - 1][A(i, 3) - 1]++;
		}

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < l; j++) {
				for (int h = 0; h < k; h++) {
					for (int f = 0; f < n; f++) {
						if (this->D3[i][h][f] < 1 || this->D33[j][h][f] < 1 || this->D4[i][j][h][f] < 1 || this->D2[h][f] < 1)
							continue;
						sum[0] += 2.0 * this->D4[i][j][h][f] * log( (this->D4[i][j][h][f] * this->D2[h][f]) / (this->D3[i][h][f] * this->D33[j][h][f]) );
                        this->D4[i][j][h][f] = 0;
					}
				}
			}
		}

		for (int i = 0; i < this->maxi; i++)
		{
			for (int j = 0; j < this->maxi; j++)
			{
				for (int k = 0; k < this->maxi; k++)
				{
                    this->D3[i][j][k] = 0;
                    this->D33[i][j][k] = 0;
				}
                this->D2[i][j] = 0;
			}
		}

		int DF = cardinality[0] * cardinality[1];
		for (int i = 2; i < cardinality.size(); i++)
			DF *= (cardinality[i] + 1);

		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;

	} else if (hDim == 5) {
		int o = cardinality[4], n = cardinality[3], k = cardinality[2], l = cardinality[1], m = cardinality[0];
//		double ***v = ThreeD(k, n, o), ****x = FourD(m, k, n, o), ****y = FourD(l, k, n, o), *****z = FiveD(m, l, k, n, o);

		for (int i = 0; i < vDim; i++) {
			this->D3[A(i, 2) - 1][A(i, 3) - 1][A(i, 4) - 1]++;
			this->D4[A(i, 0) - 1][A(i, 2) - 1][A(i, 3) - 1][A(i, 4) - 1]++;
			this->D44[A(i, 1) - 1][A(i, 2) - 1][A(i, 3) - 1][A(i, 4) - 1]++;
			this->D5[A(i, 0) - 1][A(i, 1) - 1][A(i, 2) - 1][A(i, 3) - 1][A(i, 4) - 1]++;
		}

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < l; j++) {
				for (int h = 0; h < k; h++) {
					for (int f = 0; f < n; f++) {
						for (int e = 0; e < o; e++) {
							if (this->D4[i][h][f][e] < 1 || this->D44[j][h][f][e] < 1 || this->D5[i][j][h][f][e] < 1 || this->D3[h][f][e] < 1)
								continue;
							sum[0] += 2.0 * this->D5[i][j][h][f][e] * log( (this->D5[i][j][h][f][e] * this->D3[h][f][e]) / (this->D4[i][h][f][e] * this->D44[j][h][f][e]) );
						    this->D5[i][j][h][f][e] = 0;
                        }
					}
				}
			}
		}

		for (int i = 0; i < this->maxi; i++)
		{
			for (int j = 0; j < this->maxi; j++)
			{
				for (int k = 0; k < this->maxi; k++)
				{
					for (int l = 0; l < this->maxi; l++)
					{
                        this->D4[i][j][k][l] = 0;
                        this->D44[j][j][k][l] = 0;
					}
                    this->D3[i][j][k] = 0;
				}
			}
		}
        
		int DF = cardinality[0] * cardinality[1];
		for (int i = 2; i < cardinality.size(); i++)
			DF *= (cardinality[i] + 1);

		pvalue = pchisq(sum, DF, FALSE);
		out[0] = pvalue[0];
		out[1] = sum[0];

		return out;
		
	} else {
		// int acKey, bcKey, abcKey, cKey;
		// unordered_map<int, int> Count;
		// unordered_map<int, int> ReMap;
		// int teta = 0;
		// IntegerVector tmp(hDim);

		// for (int i = 0; i < hDim; i++) {
		// 	tmp = A(i, _);
		// 	abcKey = Hash(tmp, 0, FALSE);
		// 	bcKey = Hash(tmp, 1, FALSE);
		// 	cKey = Hash(tmp, 2, FALSE);
		// 	acKey = Hash(tmp, 0, TRUE);
		// 	if (Count[abcKey] == 0) {
		// 		ReMap[teta] = abcKey;
		// 		ReMap[teta+1] = acKey;
		// 		ReMap[teta+2] = bcKey;
		// 		ReMap[teta+3] = cKey;
		// 		teta = teta + 4;
		// 		Count[abcKey]++;
		// 		Count[acKey]++;
		// 		Count[bcKey]++;
		// 		Count[cKey]++;
		// 	}
		// }

		// for (int i = 0; i < ReMap.size(); i = i + 4)
		// 	sum[0] += 2.0 * Count[ReMap[i]] == 0 * log( (Count[ReMap[i]] == 0 * Count[ReMap[i+3]] == 0) / (Count[ReMap[i+1]] == 0 * Count[ReMap[i+2]] == 0) );

		// int DF = cardinality[0] * cardinality[1];
		// for (int i = 2; i < cardinality.size(); i++)
		// 	DF *= (cardinality[i] + 1);

		// pvalue = pchisq(sum, DF, FALSE);
		// out[0] = pvalue[0];
		// out[1] = sum[0];

		return NumericVector::create(1.0, 1.0);

	}
}

IntegerMatrix MMHC::partialMatrix(const IntegerVector& pa) {
	IntegerMatrix partMat(this->vDim, pa.size());
	for (int i = 0; i < pa.size(); i++)
		partMat(_, i) = this->A(_, pa(i));
	return partMat;
}

IntegerVector MMHC::CorrespondingCardinality(const IntegerVector& pa) {
	IntegerVector tmpCardinality;
	for (int i = 0; i < pa.size(); i++)
		tmpCardinality.push_back(this->cardinality[pa[i]]);
	return tmpCardinality;
}

IntegerVector MMHC::SetCols(const IntegerVector& cpc, int T, int X) {
	IntegerVector pa;
	pa.push_back(X);
	pa.push_back(T);
	for (int i = 0; i < cpc.size(); i++)
		pa.push_back(cpc[i]);
	return pa;
}

void MMHC::UpdateCPC(List& CPC, double selected) {
	IntegerVector tmp;
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
		for (int i = 2; i < CPCLength; i++) {
			tmp = as<IntegerVector>(CPC[i]);
			tmp.push_back(selected);
			CPC.push_back(tmp);
		}
	}
}

bool MMHC::IsIn(const IntegerVector& x, double y) {
	bool out = FALSE;
	for (int i = 0; i < x.size(); i++) {
		if (x[i] == y) {
			out = TRUE;
			break;
		}
	}
	return out;
}

void MMHC::MaxMinHeuristic(const IntegerVector& cpc, List& CPCprops, IntegerVector& variablesToTest, int T) {
	IntegerMatrix statisticMatrix;
	IntegerVector pa, tmpCardinality, rejectedInLastStep = as<IntegerVector>(CPCprops[1]), temporaryMinimum = rejectedInLastStep, accepted;
	NumericVector pvalue;

	for (int i = 0; i < variablesToTest.size(); i++) {
		pa = SetCols(cpc, T, variablesToTest[i]);
		tmpCardinality = CorrespondingCardinality(pa);
		statisticMatrix = partialMatrix(pa);
		pvalue = Svalue(statisticMatrix, tmpCardinality);

		if (pvalue[0] < this->alpha) {
			if (pvalue[0] < rejectedInLastStep[1]) {
				temporaryMinimum = rejectedInLastStep;
				rejectedInLastStep = IntegerVector::create(variablesToTest[i], pvalue[0], pvalue[1]);
			} else if (pvalue[0] == rejectedInLastStep[1] && pvalue[1] > rejectedInLastStep[2]) {
				rejectedInLastStep = IntegerVector::create(variablesToTest[i], rejectedInLastStep[1], pvalue[1]);
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

void MMHC::CompatibilityToR(IntegerVector& cpc) {
	for (int i = 0; i < cpc.size(); i++)
		cpc[i]++;
}

List MMHC::Forward(int T) {
	IntegerVector cpc, variablesToTest;
	List CPC, CPCprops;
	
	IntegerVector reject = IntegerVector::create(-1.0, 1.0, 1.0), accepted = IntegerVector::create(T), tmpAccepted;
	CPCprops.push_back(T);
	CPCprops.push_back(reject);
	CPCprops.push_back(reject);

	for (int i = 0; i < this->hDim; i++)
		if (i != T)
			variablesToTest.push_back(i);

	UpdateCPC(CPC, 0);
	MaxMinHeuristic(cpc, CPCprops, variablesToTest, T);
	
	reject = as<IntegerVector>(CPCprops[1]);
	tmpAccepted = as<IntegerVector>(CPCprops[0]);

	if (reject[0] != -1.0) {
		UpdateCPC(CPC, reject[0]);
		cpc.push_back(reject[0]);
		accepted.push_back(reject[0]);
	}

	for (int i = 0; i < tmpAccepted.size(); i++)
		accepted.push_back(tmpAccepted[i]);

	variablesToTest = IntegerVector::create();
	for (int i = 0; i < this->hDim; i++)
	{
		if (!IsIn(accepted, i)) {
			variablesToTest.push_back(i);
		}
	}

	while (variablesToTest.size() > 0) {
		CPCprops[1] = IntegerVector::create(-1.0, 1.0, 1.0);
		MaxMinHeuristic(cpc, CPCprops, variablesToTest, T);

		reject = as<IntegerVector>(CPCprops[1]);
		tmpAccepted = as<IntegerVector>(CPCprops[0]);

		if (reject[0] != -1.0) {
			UpdateCPC(CPC, reject[0]);
			cpc.push_back(reject[0]);
			accepted.push_back(reject[0]);

			for (int i = 0; i < tmpAccepted.size(); i++)
				accepted.push_back(tmpAccepted[i]);

			variablesToTest = IntegerVector::create();
			for (int i = 0; i < this->hDim; i++)
			{
				if (!IsIn(accepted, i)) {
					variablesToTest.push_back(i);
				}
			}
		} else {
			variablesToTest = IntegerVector::create();
		}
	}

	return CPC;
}

IntegerVector MMHC::Backward(List& CPC, int T) {
	IntegerVector cpc = as<IntegerVector>(CPC[CPC.size()-1]), rm, pa, tmpCardinality, fromCPC;
	NumericVector pvalue;
	IntegerMatrix statisticMatrix;
	int k = -1;

	if (cpc.size() == 1) {
		return cpc;
	} else {
		for (int i = 0; i < cpc.size(); i++) {
			for (int j = 1; j < CPC.size(); j++) {
				if (j == 1) {
					fromCPC = IntegerVector::create();
				} else {
					fromCPC = as<IntegerVector>(CPC[j]);
				}

				pa = SetCols(fromCPC, T, cpc[i]);
				tmpCardinality = CorrespondingCardinality(pa);
				statisticMatrix = partialMatrix(pa);
				pvalue = Svalue(statisticMatrix, tmpCardinality);

				if (pvalue[0] > this->alpha && pvalue[0] != 1.0) {
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

void MMHC::mmpc() {
	IntegerVector cpc, pcVec, tmppc;
	List CPC, tmpPC;
	
	for (int T = 0; T < this->hDim; T++) {
		CPC = Forward(T);
		cpc = as<IntegerVector>(CPC[CPC.size()-1]);
		if (cpc.size() == 0) {
			this->pc.push_back(R_NilValue);
		} else {
			pcVec = Backward(CPC, T);
			for (int j = 0; j < pcVec.size(); j++) {
				tmpPC = Forward((int)pcVec[j]);
				tmppc = Backward(tmpPC, (int)pcVec[j]);
				if (!(IsIn(tmppc, T))) {
					for (int l = 0; l < pcVec.size(); l++) {
						if (pcVec[l] == pcVec[j]) {
							pcVec.erase(l);
						}
					}
				}
			}
			CompatibilityToR(pcVec);
			this->pc.push_back(pcVec);
		}
	}
}