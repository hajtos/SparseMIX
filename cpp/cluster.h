#pragma once

#include <set>
#include <map>
#include <list>
#include <vector>

using namespace std;

#define Point vector<int>

class Cluster {
	list< Point*> points;
	double cost;
	int* dim_counts;
	double* cost_parts;
	Point representant;
	int rep_sum;
	int dis_sum;
	int dim;
	set<int>* dim_buckets;
	double t;
	
	
	public:
	int point_count;
	Cluster(){}
	Cluster(int dim, double t, int n);
	void addPoint( Point* point);
	void removePoint( Point* point);
	double getCostWith( Point* point, bool updated=false);
	double getCostWithout( Point* point, bool updated=false);
	double getCost();
};

