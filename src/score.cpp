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

// int getN_ij();
// int getSingleN_ijk(int*, int*, int);
// double ScoreEmptyGraph(SEXP, int*, int*, double*, int*);
// double BDeu();

// int main (int argc, char *argv[]) {

// 	return 0;
// }

string int_array_to_string(NumericVector array) { //  int int_array[], int size_of_array){
	ostringstream oss("");
	for (int temp = 0; temp < array.size(); temp++)
		oss << array[temp];
	
	return oss.str();
}

// [[Rcpp::export]]
unordered_map<string, int> Unique(SEXP x) {
	NumericMatrix A(x);
	NumericVector u;
	string key;
	unordered_map<string, int> Map;

	for (int i = 0; i < A.nrow(); i++)
	{
		u = A(i, _);
		key = int_array_to_string(u);
		Map[key] = 1;
	}

	return Map;
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

int getCombiN_ijk(double *vec, double *parentVec, int dim, int j, int k) {
	int count = 0;

	for (int i = 0; i < dim; i++)
	{
		if (vec[i] == k && parentVec[i] == j)
			count++;
	}

	return count;
}

// [[Rcpp::export]]
double ScoreNodeWithNoneParent(SEXP column, SEXP N, SEXP R, SEXP Eta) {
	int *n = INTEGER(N), *r = INTEGER(R), dim = XLENGTH(column);
	double *vec = REAL(column), *eta = REAL(Eta);
	double gammaJ = *eta / 1.0, gammaK = *eta / (1.0 * (*r));
	double rScore = 0.0, qScore = 0.0;
	int n_ijk = 0, n_ij = 0;

	for (int k = 1; k <= *r; k++)
	{
		n_ijk = getSingleN_ijk(vec, dim, k);
		n_ij += n_ijk;
		rScore += lgamma( n_ijk + gammaK ) - lgamma( gammaK );
	}

	qScore += lgamma( gammaJ ) - lgamma( n_ij + gammaJ ) + rScore;


	return qScore;
}

// [[Rcpp::export]]
double ScoreNodeWithOneParent(SEXP Xi, SEXP Pa, SEXP N, SEXP R, SEXP Q, SEXP Eta) {
	int *n = INTEGER(N), *r = INTEGER(R), *q = INTEGER(Q), dim = XLENGTH(Xi);
	double *vec = REAL(Xi), *eta = REAL(Eta);
	NumericMatrix Parent(Pa);
	unordered_map<string, int> Map = Unique(Parent);
	double gammaJ = *eta / (*q), gammaK = *eta / ((*q) * (*r));
	double rScore = 0.0, qScore = 0.0;
	int n_ij, n_ijk = 0;


	for (int j = 1; j <= *q; j++)
	{
		n_ij = 0;
		rScore = 0.0;

		for (int i = 1; i <= *r; i++)
		{
			n_ij += getCombiN_ijk(vec, parentVec, dim, j, i);
		}

		for (int k = 1; k <= *r; k++)
		{
			n_ijk = getCombiN_ijk(vec, parentVec, dim, j, k);
			rScore += lgamma( n_ijk + gammaK ) - lgamma( gammaK );
		}

		qScore += lgamma( gammaJ ) - lgamma( n_ij + gammaJ ) + rScore;
	}


	return qScore;
}

// [[Rcpp::export]]
double ScoreNodeWithMultipleParent(SEXP Xi, SEXP Pa, SEXP N, SEXP R, SEXP Q, SEXP Eta, SEXP row, SEXP col) {
	int *n = INTEGER(N), *r = INTEGER(R), *q = INTEGER(Q), dim = XLENGTH(Xi);
	double *vec = REAL(Xi), *eta = REAL(Eta), *parentVec = REAL(Pa);
	double gammaJ = *eta / (*q), gammaK = *eta / ((*q) * (*r));
	double rScore = 0.0, qScore = 0.0;
	int n_ij, n_ijk = 0;

	for (int j = 1; j <= *q; j++)
	{
		n_ij = 0;
		rScore = 0.0;

		for (int i = 1; i <= *r; i++)
		{
			n_ij += getCombiN_ijk(vec, parentVec, dim, j, i);
		}

		for (int k = 1; k <= *r; k++)
		{
			n_ijk = getCombiN_ijk(vec, parentVec, dim, j, k);
			rScore += lgamma( n_ijk + gammaK ) - lgamma( gammaK );
		}

		qScore += lgamma( gammaJ ) - lgamma( n_ij + gammaJ ) + rScore;
	}


	return qScore;
}

// // [[Rcpp::export]]
// double BDeu() {
// 	unordered_map<int, int> Node;
// 	double currentScore = highestScore = 0.0;
// 	int n;

// 	for (int i = 0; i < n; i++)
// 		currentScore += (-1.0 * ScoreEmptyGraph(Mat(i, _), n, r, eta));
	
// 	while (highestScore < currentScore) {

// 		highestScore = currentScore;


// 	}
// }