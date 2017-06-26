import Cluster
import time
import random
from math import log

def nlogn(n):
    return n*log(n) if n!=0 else 0

def runAlgorithm(samples, cluster_number, dimensions, start_time=None, t=0.5, seed=0, alfa=1.0, eps=0.03):
    clusters = [Cluster.Cluster(dimensions, i, t) for i in range(cluster_number)]
    # initial cluster association, may be randomized or given
    random.seed(seed)
    #points_clusters is a list that maps each point to the cluster which it resides in
    points_clusters = [random.randint(0, cluster_number-1) for _ in samples]
    #points_clusters = [0 for _ in samples]
    for i, point in enumerate(samples):
        clusters[points_clusters[i]].addPoint(point)
    #
    updates = 1     #number of points updated in the current iteration
    iters = 0       #number of iterations performed
    energies = []
    while updates > 0:
        updates = 0
        iters += 1
        # epsilon do usuwania
        # porownanc z epsilonem(np. 3%), z alfa(0.5) i oba
        for i, point in enumerate(samples):
            mincluster = points_clusters[i]
            mincost = clusters[mincluster].cost - clusters[mincluster].getCostWithout(point)
            mincost = alfa*mincost + (1-alfa)*(nlogn(clusters[mincluster].point_count-1) \
                - nlogn(clusters[mincluster].point_count))
            for j, cluster in enumerate(clusters):
                if j == points_clusters[i] or cluster.point_count == 0:
                    continue
                diff_cost = cluster.getCostWith(point) - cluster.cost
                diff_cost = alfa*diff_cost + (1-alfa)*(nlogn(cluster.point_count)-nlogn(cluster.point_count+1))
                if diff_cost < mincost:
                    mincost = diff_cost
                    mincluster = j
            if mincluster != points_clusters[i]:
                updates += 1
                clusters[mincluster].addPoint(point)
                clusters[points_clusters[i]].removePoint(point)
                points_clusters[i] = mincluster
        for i, cl in enumerate(clusters):
            if cl.point_count < eps*len(samples) and cl.point_count != 0:
                updates += 1
                for j, point in enumerate(samples):
                    if points_clusters[j] == i:
                        z = -1
                        while True:
                            z = random.randint(0, cluster_number-1)
                            if z != i and clusters[z].point_count > 0:
                                break
                        points_clusters[j] = z
                        cl.removePoint(point)
                        clusters[z].addPoint(point)
        #print iters, updates, time.time() - start_time if start_time is not None else "", sum(cl.cost for cl in clusters)
        cost = sum([cluster.cost+nlogn(cluster.point_count) for cluster in clusters])
        energies.append(cost)
    print energies
    return (clusters, points_clusters, cost)
