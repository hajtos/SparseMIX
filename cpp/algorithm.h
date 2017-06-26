#pragma once

#include "cluster.h"
#include <vector>

void runAlgorithm( vector<Point>& samples, int cluster_number, int dimensions, int* points_clusters, double t=0.5, int seed=0, double alfa=1.0, double eps=0.03);
