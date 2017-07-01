#pragma once

#include <set>
#include <map>
#include <list>
#include <vector>

using namespace std;

#define Point vector<int>

class Cluster {
    list< Point*> points;                           //a list of the cluster's points
    double cost;                                    //the current inner cost of the cluster
    int* dim_counts;                                //the number of points containing a 1 for each position
    double* cost_parts;                             //partial values of the cost for each position
    Point representative;                           //the representative of the cluster
    int rep_sum;                                    //the number of 1s in the representative
    int dis_sum;                                    //the overall number of differences between
                                                    //the points and the representative
    int dim;                                        //the dimensionality of the data
    set<int>* dim_buckets;                          //a mapping reverse to dim_counts
    double t;                                       //the representative treshold

    public:
    int point_count;                                //the size of the cluster
    Cluster(){}                                     //default constructor, needed only for array creation
    /*
    the practical constructor:
        dim - the dimensionality of the data
        t - the representative treshold
        n - the maximal expected size(usually just the size of the dataset)
    */
    Cluster(int dim, double t, int n);
    /*
    Adds a point to the cluster, updating all neccessary values
    */
    void addPoint( Point* point);
    /*
    Removes a point from the cluster, updating all neccessary values
    */
    void removePoint( Point* point);
    /*
    Calculates the inner cost of the cluster if the given point were to be added.
    If updated is true, the cluster parameters will be updated as if the point has been added
    */
    double getCostWith( Point* point, bool updated=false);
    /*
    Calculates the inner cost of the cluster if the given point were to be removed.
    If updated is true, the cluster parameters will be updated as if the point has been removed
    */
    double getCostWithout( Point* point, bool updated=false);
    //returns the current cost of the cluster
    double getCost();
};

