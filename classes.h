#ifndef CLASSES_H
#define CLASSES_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <vector>
#include <assert.h>
#include <time.h>
using namespace std;

// vector<double> operator-(const vector<double>&,const vector<double>&);
// void operator+=(vector<double>&,const vector<double>&);
// vector<double> operator*(double,vector<double>);
// double operator|(const vector<double>&,const vector<double>&);
// double operator*(const vector<double>&,const vector<double>&);

struct MatchAndScore
{
    int match;
    int score;
};

class Matrix
{
    public:
        virtual int Size()=0;
        virtual int Get(int,int,int)=0;
        void Print();
};

class ScoreMatrix : public Matrix
{
    private:
        int dim;
        vector<vector<MatchAndScore> > A;
    public:
        ScoreMatrix();
        ScoreMatrix(int);
        ScoreMatrix(int,vector<vector<MatchAndScore> >);
        ~ScoreMatrix();

        int Size();
        int Get(int,int,int);
};

#endif