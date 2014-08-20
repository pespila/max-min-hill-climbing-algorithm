#include "mmhc.h"

unordered_map<int, int> MMHC::UniqueMap(IntegerMatrix& A) {
    IntegerVector u;
	int key;
	unordered_map<int, int> Map;
	unordered_map<int, int> Uni;
	int n = Map.size(), k = 0;

	for (int i = 0; i < this->vDim; i++)
	{
		u = A(i, _);
		key = Hash(u, 0, FALSE);
		Map[key] = 1;
		if (Map.size() != n) {
			Uni[k] = key;
			n++;
			k++;
		}
	}
	return Uni;
}

int MMHC::getSingleN_ijk(IntegerVector& vec, int k) {
	int count = 0;
	for (int i = 0; i < vec.size(); i++)
		if (vec[i] == k)
			count++;
	return count;
}

int MMHC::getVecN_ijk(IntegerVector& vec, IntegerVector& parentVec, int j, int k) {
	int count = 0;
	for (int i = 0; i < vec.size(); i++)
		if (vec[i] == k && parentVec[i] == j)
			count++;
	return count;
}

int MMHC::getMapN_ijk(IntegerVector& vec, IntegerMatrix& parentMatrix, unordered_map<int, int> parentMap, int j, int k) {
	int count = 0;
	IntegerVector u;
	for (int i = 0; i < vec.size(); i++) {
		u = parentMatrix(i, _);
		if (vec[i] == k && parentMap[j-1] == Hash(u, 0, FALSE))
			count++;
	}
	return count;
}

double MMHC::ScoreNodeWithNoneParents(IntegerVector& vec, int r) {
	double gammaJ = this->eta / 1.0, gammaK = this->eta / (1.0 * r);
	double rScore = 0.0, qScore = 0.0;
	int n_ijk = 0, n_ij = 0;

	for (int k = 1; k <= r; k++) {
		n_ijk = getSingleN_ijk(vec, k);
		n_ij += n_ijk;
		if (n_ijk != 0)
			rScore += lgamma( n_ijk + gammaK ) - lgamma( gammaK );
	}
	qScore += lgamma( gammaJ ) - lgamma( n_ij + gammaJ ) + rScore;
	return qScore;
}

double MMHC::ScoreNodeWithOneParent(IntegerVector& vec, IntegerVector& parentVec, int r, int q) {
	double gammaJ = this->eta / q, gammaK = this->eta / (q * r);
	double rScore = 0.0, qScore = 0.0;
	int n_ij, n_ijk = 0;

	for (int j = 1; j <= q; j++) {
		n_ij = 0;
		rScore = 0.0;
		for (int k = 1; k <= r; k++) {
			n_ijk = getVecN_ijk(vec, parentVec, j, k);
			n_ij += n_ijk;
			if (n_ijk != 0)
				rScore += lgamma( n_ijk + gammaK ) - lgamma( gammaK );
		}
		qScore += lgamma( gammaJ ) - lgamma( n_ij + gammaJ ) + rScore;
	}
	return qScore;
}

double MMHC::ScoreNodeWithMoreParents(IntegerVector& vec, IntegerMatrix& Parent, int r, int q) {
	unordered_map<int, int> parentMap = UniqueMap(Parent);
	double gammaJ = this->eta / q, gammaK = this->eta / ( q * r);
	double rScore = 0.0, qScore = 0.0;
	int n_ij, n_ijk;

	for (int j = 1; j <= q; j++) {
		n_ij = 0;
		rScore = 0.0;
		for (int k = 1; k <= r; k++) {
			n_ijk = getMapN_ijk(vec, Parent, parentMap, j, k);
			n_ij += n_ijk;
			if (n_ijk != 0)
				rScore += lgamma( n_ijk + gammaK ) - lgamma( gammaK );
		}
		qScore += lgamma( gammaJ ) - lgamma( n_ij + gammaJ ) + rScore;
	}
	return qScore;
}

NumericVector MMHC::ReturnParents(int i, IntegerMatrix& AdjMat) {
	NumericVector parents;
	for (int j = 0; j < this->hDim; j++)
		if (AdjMat(j, i) == 1)
			parents.push_back(j);
	return parents;
}

void MMHC::InitScore() {
	IntegerVector g;
	for (int i = 0; i < this->hDim; i++) {
		g = this->A(_, i);
		this->scores[i] = ScoreNodeWithNoneParents(g, this->cardinality[i]);
	}
}

void MMHC::ScoreGraph(IntegerMatrix& AdjMat, NumericVector& tmp) {
	IntegerMatrix parentMatrix;
	IntegerVector childVector, parentVector, allParents;
	int q;

	for (int i = 0; i < this->hDim; i++) {
		q = 1;
		allParents = ReturnParents(i, AdjMat);
		childVector = this->A(_, i);

		for (int k = 0; k < allParents.size(); k++)
			q *= this->cardinality[allParents[k]];

		if (allParents.size() == 0) {
			tmp[i] = ScoreNodeWithNoneParents(childVector, this->cardinality[i]);
		} else if (allParents.size() == 1) {
			parentVector = this->A(_, allParents[0]);
			tmp[i] = ScoreNodeWithOneParent(childVector, parentVector, this->cardinality[i], q);
		} else {
			parentMatrix = partialMatrix(allParents);
			tmp[i] = ScoreNodeWithMoreParents(childVector, parentMatrix, this->cardinality[i], q);
		}
	}
}

void MMHC::SettingEdges() {
	NumericVector tmp(this->hDim);
	double before, after;

	for (int i = 0; i < this->pc.size(); i++) {
		IntegerVector subPC = as<IntegerVector>(this->pc[i]);
		for (int j = 0; j < subPC.size(); j++) {
			if (this->graph(i, subPC[j] - 1) == 0 && this->graph(subPC[j] - 1, i) == 0) {
				this->graph(i, subPC[j] - 1) = 1;
				before = sum(this->scores);
				ScoreGraph(this->graph, tmp);
				after = sum(tmp);

				if (after > before) {
					this->scores = tmp;
				} else {
					this->graph(i, subPC[j] - 1) = 0;
				}
			}
		}
	}
}

void MMHC::AddReverseDelete(IntegerMatrix& AdjMat, NumericVector& scores) {
	srand(time(NULL));
	NumericVector tmp(this->hDim);
	double before, after;
	int rnd;

	for (int i = 0; i < this->pc.size(); i++) {
		IntegerVector subPC = as<IntegerVector>(this->pc[i]);
		for (int j = 0; j < subPC.size(); j++) {
			rnd = (rand() % 10) + 1;
			if (AdjMat(i, subPC[j] - 1) == 0 && AdjMat(subPC[j] - 1, i) == 0) {
				AdjMat(subPC[j] - 1, i) = 1;
				before = sum(scores);
				ScoreGraph(AdjMat, tmp);
				after = sum(tmp);

				if (after > before) {
					scores = tmp;
				} else {
					AdjMat(subPC[j] - 1, i) = 0;
				}
			} else if (AdjMat(i, subPC[j] - 1) == 1 && rnd > 5) {
				AdjMat(i, subPC[j] - 1) = 0;
				before = sum(scores);
				ScoreGraph(AdjMat, tmp);
				after = sum(tmp);

				if (after > before) {
					scores = tmp;
				} else {
					AdjMat(i, subPC[j] - 1) = 1;
				}
			}  else if (AdjMat(subPC[j] - 1, i) == 1 && rnd <= 5) {
				AdjMat(subPC[j] - 1, i) = 0;
				AdjMat(i, subPC[j] - 1) = 1;
				before = sum(scores);
				ScoreGraph(AdjMat, tmp);
				after = sum(tmp);

				if (after > before) {
					scores = tmp;
				} else {
					AdjMat(subPC[j] - 1, i) = 1;
					AdjMat(i, subPC[j] - 1) = 0;
				}
			}
		}
	}
}


void MMHC::mmhc() {
	NumericVector tmpScores;
	IntegerMatrix tmpAdjMat(this->hDim, this->hDim);

	InitScore();
	SettingEdges();
	AddReverseDelete(tmpAdjMat, tmpScores);

	int count = 1;
	for (int i = 0; i < 1000; i++) {
		if (count == 15)
			break;

		tmpScores = this->scores;
		tmpAdjMat = this->graph;
		AddReverseDelete(tmpAdjMat, tmpScores);
		
		if (sum(this->scores) > sum(tmpScores)) {
			count++;
		} else {
			this->graph = tmpAdjMat;
			this->scores = tmpScores;
			count = 1;
		}
	}
	this->score = sum(this->scores);
}