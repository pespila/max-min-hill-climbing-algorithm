#include <iostream>
#include <omp.h>
#include <unordered_map>
#include <vector>
#include <time.h>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

string int_array_to_string(vector<int> array) { //  int int_array[], int size_of_array){
	ostringstream oss("");
	for (int temp = 0; temp < array.size(); temp++)
		oss << array[temp];

	return oss.str();
}

int main(void) {
	srand (time(NULL));
	vector<int> x(10);

	for (int i = 0; i < x.size(); i++)
	{
		x[i] = rand()%3;
		cout << x[i] << endl;
	}

	cout << endl << endl;

	string out;

	out = int_array_to_string(x);
	

	vector<int> y(10);
	for (int i = 0; i < y.size(); i++)
	{
		y[i] = out[i] - '0';
		cout << y[i] << endl;
	}
	
	unordered_map<string, int> Map;
	Map[out] = 6;
	return 0;
}