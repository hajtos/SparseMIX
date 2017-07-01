#include "algorithm.h"
#include "cluster.h"
#include <iostream>
#include <cmath>
#include <stdlib.h>

//utility function
static double nlogn(int n) {
    return n==0? 0 : log(n)*n;
}

void runAlgorithm( vector<Point>& samples, int cluster_number, int dimensions,
        int* points_clusters, double t, int seed, double alfa, double eps) {
    Cluster* clusters = new Cluster[cluster_number];
    int i, mincluster;
    int updates = 1;
    int iters = 0;
    double mincost, diff_cost;
    srand(seed);

    for (i=0;i<cluster_number;++i)
        clusters[i] = Cluster(dimensions, t, samples.size());
    // initial cluster association
    for (i=0;i<samples.size();++i) {
        points_clusters[i] = rand()%cluster_number;
        clusters[points_clusters[i]].addPoint(&samples[i]);
    }
    while (updates > 0) { //main loop
        updates = 0;
        iters += 1;
        for (i=0; i<samples.size();++i) { //try reassigning each point
            mincluster = points_clusters[i];
            Point* point = &samples[i];
            //calculate the difference between the cost of the
            //current point's cluster and of that point were to be removed
            mincost = clusters[mincluster].getCost() -
                clusters[mincluster].getCostWithout(point);
            mincost = alfa*mincost + (1-alfa)*(nlogn(clusters[mincluster].point_count-1)
                - nlogn(clusters[mincluster].point_count));
            for (int j=0; j<cluster_number;++j) {
                Cluster& cluster = clusters[j];
                if (j == points_clusters[i])
                    continue;
                //calculate an analogous difference for each cluster and find the smallest
                diff_cost = cluster.getCostWith(point) - cluster.getCost();
                diff_cost = alfa*diff_cost + (1-alfa)*(nlogn(cluster.point_count)
                    -nlogn(cluster.point_count+1));
                if (diff_cost < mincost) {
                    mincost = diff_cost;
                    mincluster = j;
                }
            }
            if (mincluster != points_clusters[i]) {
                //if the best cluster is other than the previous, reassign
                updates += 1;
                clusters[mincluster].addPoint(point);
                clusters[points_clusters[i]].removePoint(point);
                points_clusters[i] = mincluster;
            }
        }
        for (i=0; i<cluster_number;++i) { //reduction of too small clusters
            if (clusters[i].point_count < eps*samples.size()
                    && clusters[i].point_count != 0) {
                updates += 1;
                for (int j=0; j<samples.size();++j) {
                    //reassign each point to a random non-empty cluster
                    Point* point = &samples[j];
                    if (points_clusters[j] == i) {
                        int z = i;
                        while (1) {
                            z = rand()%cluster_number;
                            if (z != i && clusters[z].point_count > 0)
                                break;
                        }
                        points_clusters[j] = z;
                        clusters[i].removePoint(point);
                        clusters[z].addPoint(point);
                    }
                }
            }
        } //end small clusters reduction
    } //end main loop
    delete[] clusters;
}
