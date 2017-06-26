// codeclusterer.cpp : Defines the entry point for the console application.
//



#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
using namespace std;

#include "algorithm.h"
#include "cluster.h"


int main(int argc, char* argv[]) {
	if (argc < 4) {
		cout << "Usage: ./main <filename> <no. of clusters> <T parameter> [<seed>]" << endl;
		return 1;
	}

	string filename = argv[1];
	int cl_num = atoi(argv[2]);
	float T = atof(argv[3]);
	int seed = 0;
	if (argc > 4) {
		seed = atoi(4);
	}
	ifstream input(filename.c_str());

	int size, dimensions;
	input >> size >> dimensions;

	string buffer;
	getline(input, buffer);
	vector<Point> samples;

	for (int i = 0; i<size; ++i) {
		getline(input, buffer);
		istringstream ss(buffer);

		int val;
		Point point;
		while (ss.good()) {
			ss >> val;
			point.push_back(val);
		}
		samples.push_back(point);
	}
	int* points_clusters = new int[size];

	runAlgorithm(samples, cl_num, dimensions, points_clusters, T, seed, 0.5, 0.03);

	for (int i = 0; i<size; ++i) {
		cout << points_clusters[i] << " ";
	}
	cout << "\b" << endl;

	return 0;
}



