#pragma once

#include "cluster.h"
#include <vector>

/*
    performs the SparseMIX algorithm

    samples - a vector of points to be clustered, each point being a vector of positions of 1s(integers)
    cluster_number - the number of clusters to partition data into
    dimensions - the number of dimensions in the data
    point_clusters - output parameter, an array that will contain the cluster numbers assigned to each point,
        needs to be allocated before function call
    t - the threshold of the method after which a position is added to the representative
    seed - the seed for the randomization
    alfa - the percentage of the cost function attributed to the inner cluster cost. The rest of the cost
        will be composed of the cluster coding (n*log(n) where n is the cluter size). A value of 1 disables
        cluster coding in the cost completely.
    eps - if any cluster is reduced to less than eps*samples.size() points, it will be reduced completely
*/
void runAlgorithm( vector<Point>& samples, int cluster_number, int dimensions, 
    int* points_clusters, double t=0.5, int seed=0, double alfa=1.0, double eps=0.03);
