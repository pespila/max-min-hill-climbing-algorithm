#include "classes.h"

void Matrix::Print() {
    if(Size()<=10) {
        for (int i = 0; i < Size(); i++)
        {
            for (int j = 0; j < Size(); j++)
            {
                printf("%d ", Get(i,j,0));
            }
            printf("\n");
        }
        printf("\n");

        for (int i = 0; i < Size(); i++)
        {
            for (int j = 0; j < Size(); j++)
            {
                printf("%d ", Get(i,j,1));
            }
            printf("\n");
        }
        printf("\n");
    }
}

ScoreMatrix::ScoreMatrix(){}

ScoreMatrix::ScoreMatrix(int dim){
    srand(time(NULL));
    this->dim=dim;
    struct MatchAndScore MAS;
    vector<MatchAndScore> tmp(this->dim,MAS);
    A.assign(this->dim,tmp);

    for (int i = 0; i < Size(); i++)
    {
        for (int j = 0; j < Size(); j++)
        {
            if(i==j) {
                A[i][j].match=0;
                A[i][j].score=0;

            }
            else {
                A[i][j].match=rand()%2;
                A[i][j].score=rand()%10;
            }
        }
    }
}

ScoreMatrix::~ScoreMatrix() {
    vector<vector<MatchAndScore> >().swap(A);
}

ScoreMatrix::ScoreMatrix(int dim,vector<vector<MatchAndScore> > A) {
    this->dim=dim;
    this->A.resize(this->dim);
    this->A=A;
}

int ScoreMatrix::Get(int i,int j,int k) {
    if(k==0) return this->A[i][j].match;
    else return this->A[i][j].score;
}

int ScoreMatrix::Size() {
    return this->dim;
}