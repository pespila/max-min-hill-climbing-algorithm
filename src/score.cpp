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

// [[Rcpp::export]]
int getR(NumericVector& f) {
	unordered_map<int, int> Map;

	for (NumericVector::iterator it = f.begin(); it != f.end(); it++) {
		Map[*it] = 1;
	}

	return Map.size();
}

string int_array_to_string(NumericVector& array) { //  int int_array[], int size_of_array){
	ostringstream oss("");
	for (int temp = 0; temp < array.size(); temp++)
		oss << array[temp];
	
	return oss.str();
}

// [[Rcpp::export]]
unordered_map<int, string> UniqueMap(SEXP x) {
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
		Map[key] = 1;
		if (Map.size() != n) {
			Uni[k] = key;
			n++;
			k++;
		}
	}

	return Uni;
}

int getSingleN_ijk(double *vec, int dim, int k) {
	int count = 0;

	for (int i = 0; i < dim; i++)
	{
		if (vec[i] == k)
			count++;
	}

	return count;
}

int getVecN_ijk(double *vec, double *parentVec, int dim, int j, int k) {
	int count = 0;

	for (int i = 0; i < dim; i++)
	{
		if (vec[i] == k && parentVec[i] == j)
			count++;
	}

	return count;
}

int getMapN_ijk(double *vec, NumericMatrix& parentMatrix, unordered_map<int, string> parentMap, int dim, int j, int k) {
	int count = 0;
	NumericVector u;

	for (int i = 0; i < dim; i++)
	{
		u = parentMatrix(i, _);
		if (vec[i] == k && parentMap[j-1] == int_array_to_string(u))
			count++;
	}

	return count;
}

// [[Rcpp::export]]
double ScoreNodeWithNoneParents(SEXP column, SEXP N, int r, double eta) {
	int *n = INTEGER(N), dim = XLENGTH(column);
	double *vec = REAL(column);
	double gammaJ = eta / 1.0, gammaK = eta / (1.0 * r);
	double rScore = 0.0, qScore = 0.0;
	int n_ijk = 0, n_ij = 0;

	for (int k = 1; k <= r; k++)
	{
		n_ijk = getSingleN_ijk(vec, dim, k);
		n_ij += n_ijk;
		if (n_ijk != 0)
			rScore += lgamma( n_ijk + gammaK ) - lgamma( gammaK );
	}

	qScore += lgamma( gammaJ ) - lgamma( n_ij + gammaJ ) + rScore;


	return qScore;
}

// [[Rcpp::export]]
double ScoreNodeWithOneParent(SEXP Xi, SEXP Pa, SEXP N, int r, int q, double eta) {
	int *n = INTEGER(N), dim = XLENGTH(Xi);
	double *vec = REAL(Xi), *parentVec = REAL(Pa);
	double gammaJ = eta / q, gammaK = eta / (q * r);
	double rScore = 0.0, qScore = 0.0;
	int n_ij, n_ijk = 0;

	for (int j = 1; j <= q; j++)
	{
		n_ij = 0;
		rScore = 0.0;

		for (int k = 1; k <= r; k++)
		{
			n_ijk = getVecN_ijk(vec, parentVec, dim, j, k);
			n_ij += n_ijk;
			if (n_ijk != 0)
				rScore += lgamma( n_ijk + gammaK ) - lgamma( gammaK );
		}

		qScore += lgamma( gammaJ ) - lgamma( n_ij + gammaJ ) + rScore;
	}


	return qScore;
}

// [[Rcpp::export]]
double ScoreNodeWithMoreParents(SEXP Xi, SEXP Pa, SEXP N, int r, int q, double eta) {
	int *n = INTEGER(N), dim = XLENGTH(Xi);
	double *vec = REAL(Xi);
	NumericMatrix Parent(Pa);
	unordered_map<int, string> parentMap = UniqueMap(Parent);
	double gammaJ = eta / q, gammaK = eta / ( q * r);
	double rScore = 0.0, qScore = 0.0;
	int n_ij, n_ijk;
	NumericVector N_ijk(r);

	for (int j = 1; j <= q; j++)
	{
		n_ij = 0;
		rScore = 0.0;

		for (int k = 1; k <= r; k++)
		{
			n_ijk = getMapN_ijk(vec, Parent, parentMap, dim, j, k);
			n_ij += n_ijk;
			if (n_ijk != 0)
				rScore += lgamma( n_ijk + gammaK ) - lgamma( gammaK );
		}

		qScore += lgamma( gammaJ ) - lgamma( n_ij + gammaJ ) + rScore;
	}


	return qScore;
}

NumericVector ReturnParents(int i, NumericMatrix& AdjMat, SEXP N) {
	NumericVector parents;
	int *n = INTEGER(N);
	for (int j = 0; j < *n; j++)
	{
		if (AdjMat(j, i) == 1) {
			parents.push_back(j);
		}
	}
	return parents;
}

NumericVector InitScore(NumericMatrix& A, SEXP N, NumericVector& R, double eta) {
	int *n = INTEGER(N);
	NumericVector g, score(*n);

	for (int i = 0; i < *n; i++) {
		g = A(_, i);
		score[i] = ScoreNodeWithNoneParents(g, N, R[i], eta);
	}

	return score;
}

// [[Rcpp::export]]
NumericMatrix partialMatrix(const NumericMatrix& A, NumericVector& pa) {
	NumericMatrix partMat(A.nrow(), pa.size());

	for (int i = 0; i < pa.size(); i++)
	{
		partMat(_, i) = A(_, pa(i));
	}

	return partMat;
}

NumericVector ScoreGraph(const NumericMatrix& A, NumericMatrix& AdjMat, SEXP N, const NumericVector& R, double eta) {
	NumericVector childVector, parentVector, allParents, score;
	NumericMatrix parentMatrix;
	int q;

	for (int i = 0; i < AdjMat.ncol(); i++)
	{
		q = 1;
		allParents = ReturnParents(i, AdjMat, N);
		childVector = A(_, i);

		for (int k = 0; k < allParents.size(); k++)
			q *= R[allParents[k]];

		if (allParents.size() == 0) {
			score.push_back(ScoreNodeWithNoneParents(childVector, N, R[i], eta));
		} else if (allParents.size() == 1) {
			parentVector = A(_, allParents[0]);
			score.push_back(ScoreNodeWithOneParent(childVector, parentVector, N, R[i], q, eta));
		} else {
			parentMatrix = partialMatrix(A, allParents);
			score.push_back(ScoreNodeWithMoreParents(childVector, parentMatrix, N, R[i], q, eta));
		}
	}

	return score;
}

void SettingEdges(const NumericMatrix& A, NumericMatrix& AdjMat, NumericVector& scores, const List& PC, SEXP N, const NumericVector& R, double eta) {
	NumericVector childVector, parentVector, allParents, tmp;
	NumericMatrix parentMatrix;
	double add, before, after;
	int q;

	for (int i = 0; i < PC.size(); i++)
	{
		NumericVector pc = as<NumericVector>(PC[i]);
		// q = 1;
		// add = scores[i];

		for (int j = 0; j < pc.size(); j++) {

			if (AdjMat(i, pc[j] - 1) == 0 && AdjMat(pc[j] - 1, i) == 0) {
				
				AdjMat(i, pc[j] - 1) = 1;

				before = sum(scores);
				tmp = ScoreGraph(A, AdjMat, N, R, eta);
				after = sum(tmp);

				if (after > before) {
					scores = tmp;
				} else {
					AdjMat(i, pc[j] - 1) = 0;
				}
			}
		}
	}
}

void AddReverseDelete(const NumericMatrix& A, NumericMatrix& AdjMat, NumericVector& scores, const List& PC, SEXP N, const NumericVector& R, double eta) {
	srand(time(NULL));
	NumericVector tmp;
	double before, after;
	int rnd;

	for (int i = 0; i < PC.size(); i++)
	{
		NumericVector pc = as<NumericVector>(PC[i]);

		for (int j = 0; j < pc.size(); j++) {

			rnd = (rand() % 10) + 1;

			if (AdjMat(i, pc[j] - 1) == 0 && AdjMat(pc[j] - 1, i) == 0) {
				
				AdjMat(pc[j] - 1, i) = 1;

				before = sum(scores);
				tmp = ScoreGraph(A, AdjMat, N, R, eta);
				after = sum(tmp);

				if (after > before) {
					scores = tmp;
				} else {
					AdjMat(pc[j] - 1, i) = 0;
				}
			} else if (AdjMat(i, pc[j] - 1) == 1 && rnd > 7) {
				
				AdjMat(i, pc[j] - 1) = 0;

				before = sum(scores);
				tmp = ScoreGraph(A, AdjMat, N, R, eta);
				after = sum(tmp);

				if (after > before) {
					scores = tmp;
				} else {
					AdjMat(i, pc[j] - 1) = 1;
				}
			}  else if (AdjMat(pc[j] - 1, i) == 1 && rnd <= 7) {

				AdjMat(pc[j] - 1, i) = 0;
				AdjMat(i, pc[j] - 1) = 1;

				before = sum(scores);
				tmp = ScoreGraph(A, AdjMat, N, R, eta);
				after = sum(tmp);

				if (after > before) {
					scores = tmp;
				} else {
					AdjMat(pc[j] - 1, i) = 1;
					AdjMat(i, pc[j] - 1) = 0;
				}
			}
		}
	}
}

// [[Rcpp::export]]
SEXP BDeu(SEXP x, SEXP y, SEXP z, SEXP N) {
	int *n= INTEGER(N);
	double eta = 2.2;
	NumericVector R(y), scores, tmpScores;
	NumericMatrix A(x), AdjMat(*n, *n), tmpAdjMat(*n, *n);
	List PC(z);

	scores = InitScore(A, N, R, eta);
	SettingEdges(A, AdjMat, scores, PC, N, R, eta);
	AddReverseDelete(A, tmpAdjMat, tmpScores, PC, N, R, eta);
	
	for (int i = 0; i < 5; i++)
	{
		tmpScores = scores;
		tmpAdjMat = AdjMat;
		AddReverseDelete(A, tmpAdjMat, tmpScores, PC, N, R, eta);
		
		if (sum(scores) >= sum(tmpScores)) {
			break;
		} else {
			AdjMat = tmpAdjMat;
			scores = tmpScores;
		}
	}

	// for (int i = 0; i < 100; i++)
	// {
	// 	if (count == 2)
	// 		break;

	// 	tmpScores = scores;
	// 	tmpAdjMat = AdjMat;
	// 	AddReverseDelete(A, tmpAdjMat, tmpScores, PC, N, R, eta);
		
	// 	if (sum(scores) <= sum(tmpScores)) {
	// 		AdjMat = tmpAdjMat;
	// 		scores = tmpScores;
	// 		count++;
	// 	} else {
	// 		count = 1;
	// 	}
	// }

	cout << sum(scores) << endl;
	
	return AdjMat;
}