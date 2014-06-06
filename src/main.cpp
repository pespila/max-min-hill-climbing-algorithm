#include <iostream>
#include <omp.h>
// #include <unordered_map>
#include <vector>
using namespace std;

int main(void) {
	#pragma omp parallel
	cout << omp_get_thread_num() << endl;

	return 0;
}

// int main(void) {
// 	typedef unordered_map<int*,int> map_t;
// 	map_t Map;
// 	// vector<double> x(5), y(5, 1.0), z(5, 2.0);
// 	int *x = (int*)malloc(5*sizeof(int));
// 	int *y = (int*)malloc(5*sizeof(int));
// 	int *z = (int*)malloc(5*sizeof(int));
// 	int *tmp = (int*)malloc(5*sizeof(int));

// 	vector<vector<int> > A;
// 	vector<int> temp(100);
// 	A.assign(5, temp);

// 	int n = 1, m = 2;

// 	for (int i = 0; i < 5; i++)
// 	{
// 		for (int j = 0; j < 100; j++)
// 		{
// 			if((i+j)%2 == 0)
// 				A[i][j] = n;
// 			else {
// 				A[i][j] = m;
// 			}
// 		}
// 	}

// 	for (int i = 0; i < 100; i++)
// 	{
// 		for (int j = 0; j < 5; j++)
// 		{
// 			tmp[j] = A[j][i];
// 		}
// 		Map[tmp] = 1;
// 		// Map.insert(map_t::value_type(tmp, 1));
// 	}

// 	cout << Map.size() << endl;

// 	return 0;
// }