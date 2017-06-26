

#include "algorithm.h"
#include "cluster.h"
#include <iostream>
#include <cmath>

static double nlogn(int n) {
	return n==0? 0 : log(n)*n;
}

void runAlgorithm( vector<Point>& samples, int cluster_number, int dimensions, int* points_clusters, double t, int seed, double alfa, double eps) {
	Cluster* clusters = new Cluster[cluster_number];
	int i, mincluster;
	int updates = 1;
	int iters = 0;
	double mincost, diff_cost;

	for (i=0;i<cluster_number;++i)
		clusters[i] = Cluster(dimensions, t, samples.size());
    cout.flush();
    // initial cluster association, may be randomized or given
    //random.seed(seed)
    //points_clusters is a list that maps each point to the cluster which it resides in
    //points_clusters = [random.randint(0, cluster_number-1) for _ in samples]
    for (i=0;i<samples.size();++i) {
    	points_clusters[i] = i%cluster_number;
    	clusters[points_clusters[i]].addPoint(&samples[i]);
    	double d = clusters[points_clusters[i]].getCost();
	}
    while (updates > 0) {
        updates = 0;
        iters += 1;
        for (i=0; i<samples.size();++i) {
            mincluster = points_clusters[i];
             Point* point = &samples[i];
            mincost = clusters[mincluster].getCost() - clusters[mincluster].getCostWithout(point);
            mincost = alfa*mincost + (1-alfa)*(nlogn(clusters[mincluster].point_count-1)
                - nlogn(clusters[mincluster].point_count));
            for (int j=0; j<cluster_number;++j) {
            	Cluster& cluster = clusters[j];
                if (j == points_clusters[i])
                    continue;
                diff_cost = cluster.getCostWith(point) - cluster.getCost();
                diff_cost = alfa*diff_cost + (1-alfa)*(nlogn(cluster.point_count)-nlogn(cluster.point_count+1));
                if (diff_cost < mincost) {
                    mincost = diff_cost;
                    mincluster = j;
				}
			}
            if (mincluster != points_clusters[i]) {
                updates += 1;
                clusters[mincluster].addPoint(point);
                clusters[points_clusters[i]].removePoint(point);
                points_clusters[i] = mincluster;
            }
		}
        for (i=0; i<cluster_number;++i) {
            if (clusters[i].point_count < eps*samples.size() && clusters[i].point_count != 0) {
                updates += 1;
        		for (int j=0; j<samples.size();++j) {
        			 Point* point = &samples[j];
                    if (points_clusters[j] == i) {
                        int z = i;
                        while (1) {
                            z = z+1;
                            if (z != i && clusters[z].point_count > 0)
                                break;
						}
                        points_clusters[j] = z;
                        clusters[i].removePoint(point);
                        clusters[z].addPoint(point);
					}
				}
			}
		}
	}
    delete[] clusters;

    //cost = sum([cluster.cost+nlogn(cluster.point_count) for cluster in clusters])
    //return (clusters, points_clusters, cost)*/
}
