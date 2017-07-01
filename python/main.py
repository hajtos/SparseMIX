import Algorithm
import sys

def main(argv):
    if len(argv) < 4:
        print "Usage: python main.py <filename> <no. of clusters> <T parameter> [<seed>]"
    filename = argv[1]
    T = argv[2]
    clust_num = argv[3]
    seed = 0
    if len(argv)>4:
        seed = argv[4]
    data = open(filename).read().split("\n")[:-1]
    dimension = int(data[0].split(" ")[1])
    data = [set([int(v) for v in dat.split(" ") if v.strip()]) for dat in data[1:]]
    clusters, points_clusters, cost = Algorithm.runAlgorithm(data, clust_num, dimension, alfa=0.5, seed=seed, t=T)
    print " ".join(point_clusters)

if __name__ == "__main__":
    main(sys.argv)
