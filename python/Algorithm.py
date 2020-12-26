from typing import List, Set, Tuple

from python import Cluster
import random
from math import log

def nlogn(n: float) -> float:
    return n*log(n) if n!=0 else 0

def runAlgorithm(samples: List[Set[int]], cluster_number: int, dimensions: int, t: float=0.5, seed: int=0,
                 alfa: float=1.0, eps: float=0.03) -> Tuple[List[Cluster.Cluster], List[int], float]:
    """
    performs the SparseMIX algorithm

    samples - a list of points to be clustered, each point being a set of non-zero positions(integers)
    cluster_number - the number of clusters to partition data into
    dimensions - the number of dimensions in the data
    t - the threshold of the method after which a position is added to the representative
    seed - the seed for the randomization
    alfa - the percentage of the cost function attributed to the inner cluster cost. The rest of the cost
        will be composed of the cluster coding (n*log(n) where n is the cluter size). A value of 1 disables
        cluster coding in the cost completely.
    eps - if any cluster is reduced to less than eps*samples.size() points, it will be reduced completely
    """
    clusters = [Cluster.Cluster(dimensions, i, t) for i in range(cluster_number)]
    # initial cluster association
    random.seed(seed)
    points_clusters = [random.randint(0, cluster_number-1) for _ in samples]
    for i, point in enumerate(samples):
        clusters[points_clusters[i]].addPoint(point)
    updates = 1     #number of points updated in the current iteration
    iters = 0       #number of iterations performed
    while updates > 0: #main loop
        updates = 0
        iters += 1
        for i, point in enumerate(samples): #try reassigning each point
            mincluster = points_clusters[i]
            #calculate the difference between the cost of the
            #current point's cluster and of that point were to be removed
            mincost = clusters[mincluster].cost - clusters[mincluster].getCostWithout(point)
            mincost = alfa*mincost + (1-alfa)*(nlogn(clusters[mincluster].point_count-1) \
                - nlogn(clusters[mincluster].point_count))
            for j, cluster in enumerate(clusters):
                if j == points_clusters[i] or cluster.point_count == 0:
                    continue
                #calculate an analogous difference for each cluster and find the smallest
                diff_cost = cluster.getCostWith(point) - cluster.cost
                diff_cost = alfa*diff_cost + (1-alfa)*(nlogn(cluster.point_count)-nlogn(cluster.point_count+1))
                if diff_cost < mincost:
                    mincost = diff_cost
                    mincluster = j
            if mincluster != points_clusters[i]:
                #if the best cluster is other than the previous, reassign
                updates += 1
                clusters[mincluster].addPoint(point)
                clusters[points_clusters[i]].removePoint(point)
                points_clusters[i] = mincluster
        for i, cl in enumerate(clusters): #reduction of too small clusters
            if cl.point_count < eps*len(samples) and cl.point_count != 0:
                updates += 1
                for j, point in enumerate(samples):
                    #reassign each point to a random non-empty cluster
                    if points_clusters[j] == i:
                        z = -1
                        while True:
                            z = random.randint(0, cluster_number-1)
                            if z != i and clusters[z].point_count > 0:
                                break
                        points_clusters[j] = z
                        cl.removePoint(point)
                        clusters[z].addPoint(point)
        cost = sum([cluster.cost+nlogn(cluster.point_count) for cluster in clusters])
    return (clusters, points_clusters, cost)
